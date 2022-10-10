#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


Main file to put together a coherent dataset for RL modeling of hypotension management
Code adapted from joseph futoma (https://github.com/dtak/POPCORN-POMDP) by anonymous author of relative sparsity manuscript
Processes csv files created by bash scripts in the above github directory
Note that anonymous rolled back to fixed interval decision points, focused only on MICU,  and
that trajectory starts at first hypotensive measurement instead of at beginning of stay.


"""

import numpy as np
import pandas as pd
import os
import sys
from datetime import datetime
from datetime import timedelta
import pickle
from time import time
import pdb


# use urine output not just bp in reward
urine_in_reward=0

# use % change in map as reward
delta_rew = 1

# use last bp as reward
last_map_rew = 0

def get_first_hypo(this_maps):
	ltSICKTHRESH=this_maps[this_maps.valuenum<=MAP_SICK_THRESH]

	if ltSICKTHRESH.empty:
		print("no hyp data")
		pdb.set_trace()
		#first_low_ix=0
		#first_low_time_line = this_maps.iloc[[0]]
		#start_time = first_low_time_line.charttime.values[0]

	else:
		first_low_ix = ltSICKTHRESH.index[0]
		first_line =  ltSICKTHRESH.iloc[[0]]
		first_low_time = first_line.charttime
		# not sure why this start time value needs 
		# to be taken like this, but else won't broadcast
		start_time = first_low_time.values[0]
	return (start_time,first_low_ix)

#comp_tag='/Users/anonymous/'
print("make sure path correct for computer")

ctag='/Users/anonymous'
loc = '/Box/MIMIC/mimic-iii-v1-4/'
sys.path.insert(1, ctag+loc+'extract-scripts')

from data_load_utils import load_labs_and_vitals

#when you play around with data locally set this to path where everything is stored
#PATH_TO_QUERY_DATA = "/n/dtak/mimic-pipeline-tools/mimic-iii-v1-4/query-data/"
PATH_TO_QUERY_DATA = ctag+loc+"query-data2/"
path_to_query_data = ctag+loc+"query-data2/"

#PATH_TO_CLEANED_DATA = "/n/dtak/mimic-pipeline-tools/hypotension/model-data/"
PATH_TO_CLEANED_DATA = ctag+loc+"hypotension-RL/model-data2/"
path_to_model_data= ctag+loc+"hypotension-RL/model-data2/"

pd.set_option("display.max_columns",101)
pd.set_option("display.max_rows",101)

###########
#params for dataset construction

#default thresholds we'll use throughout, so not explicitly in filenames (for brevity)
# Note, this threshold is also in mimic-iii-v1-4/extract-scripts/get_cohort_baseline_info.py 
LOS_LOWER_THRESH = 0#300#12#12  #exclude if LOS is shorter than this many hours


# to decide when to start hypotension
MAP_SICK_THRESH = 60 #criteria for "sick"...for now... # anonymous added use for cohort select, time disc, and for start of traj


#this is a much harsher filter, may want to play around with this one...
# this filter is only applied to files with 'Xbpleq65' in the filename
# anonymous: is this used? This is not present in this file. so I think no
# NOT USED
MAP_FILTER_THRESH = 65
MAP_NUM_BELOW_THRESH = 1 #must have at least this many MAP values below MAP_FILTER_THRESH to be retained


# There will be (LOS_CAP/TIME_GRID_INT_HRS - 1) decision point
LOS_CAP = 20#48 #last time point considered, if longer than this cap here & end traj. note this is in hours

# anonymous: Joe had used a decision point method based on when people took decisions, but 
# we returned  to grid. Grid will be according to these intervals
# I guess IV norepi has peak effect  at like 2.5 min


TIME_GRID_INT_HRS = 1/4

#TODO: eventually later on might implement different args for some of these design choices
BP_AGG_FUNC = 'min' #'min'; how to aggregate BP values to get a single value...
BASELINE_IMPUTE_METHOD = 'median' #for now the only option is population median, but maybe improve later
TIMESERIES_IMPUTE_METHOD = 'last' #for now the only option is sample-and-hold or last-one-carry-forward

# start w first hypo?
hypo_start=1


SEED = 8675309 #seed for reproducing any results that involve randomness
np.random.seed(SEED)

##########
########## load in the baseline/static data for the cohort
##########

cohort_dat = pd.read_csv(PATH_TO_QUERY_DATA+'cohort.csv')
cohort_dat = cohort_dat.sort_values('ICUSTAY_ID')
# do any additional filtering beyond cohort construction, if desired...

#anonymous: can maybe filter more here.

cohort_dat = cohort_dat.loc[cohort_dat['LOS'] >= LOS_LOWER_THRESH,:]

# anonymous added this, since final was defined later.
# It is still filtered later after all_raw_action data and all_state_data are created
# The second filter I moved to after creation of all_state_data, the first filter I added
# since later in the code there was a keyerror for one of the ids not in all_rawaction data.
# I think he is just cutting out ones w few measurements or actions

# anonymous:restrict to MICU
cohort_dat = cohort_dat.loc[cohort_dat.FIRST_CAREUNIT=='MICU',]
#cohort_dat=cohort_dat.iloc[1:1000]
nrows_to_read = None#10000

all_ICU_IDs = np.array(cohort_dat['ICUSTAY_ID'])


# anonymous: if they are never hypotensive, don't keep them.  This will
# also remove anyone with no map measurements. those people not hypo anyway
# We can do this because we are going to start our trajectory at the first hypotensive measurement
# so our stategy starts at hypotension observation


##########
########## LOAD all time series data
##########

all_ts_dats, ts_medians, all_ts_vars = load_labs_and_vitals(
	PATH_TO_QUERY_DATA,cohort_dat,nrows_to_read)

print('time series variables all loaded')


# map dat is sorted by icu id and chartime
map_dat = all_ts_dats['map']#pd.read_csv(PATH_TO_QUERY_DATA+'map.csv')

# exclude patients who had no hypotension.
# we can do this because we will start at hypotension onset
# maybe further exclude patients who don t have enough measure
# after hyp onset?  We are further going to start the problem
# directly at hypotension onset, ignoring measurements before
# we will lose then info on whether they were previously given 
# vaso

hyp_IDs = []
for ID in all_ICU_IDs:
	this_map_dat = map_dat.loc[map_dat.icustay_id==ID,]
	if (this_map_dat.valuenum<MAP_SICK_THRESH).sum()>0:
		hyp_IDs.append(ID)

print("hyp",len(hyp_IDs))
# possibly cut trajectories to start after hypotension.

# remove patients without hypo from map 
#map_dat =  map_dat.loc[map_dat.icustay_id.isin(hyp_IDs),]
#start_times = {}
meas_after_hyps = []
for ID in hyp_IDs:
 	this_map_dat = map_dat.loc[map_dat.icustay_id==ID,]
 	start_time,_ = get_first_hypo(this_map_dat)
 	#start_times[ID] = start_time
 	after_hyp = this_map_dat.charttime>=start_time
 	this_maps_after_hyp = this_map_dat.loc[after_hyp,]
 	if len(this_maps_after_hyp)>1:
 		meas_after_hyps.append(ID)
print("hyp and meas",len(meas_after_hyps))

# 	after_hyps.append(after_hyp)
# after_hyps = [item for sublist in after_hyps for item in sublist]

# map_dat = map_dat.loc[after_hyps,]
# # only take measurements after hyp, for each stay

# # I think it's ok to just do this for states. those
# # are the only ones agg over time. Fluids and vaso I guess we will still get
# # namely, we do not then take all the maps before hypo onset
# pdb.set_trace()
# for k in  all_ts_dats.keys():
# 	all_ts_dats[k] = all_ts_dats[k].loc[after_hyps,]
# 	#


# make sure have measurements after hypo

cohort_dat =  cohort_dat.loc[cohort_dat.ICUSTAY_ID.isin(meas_after_hyps),]
final_ICU_IDs = cohort_dat.ICUSTAY_ID

# anonymous added this, pop_dat wasn't defined, but it is used much below.


pop_dat = cohort_dat.copy()
pop_dat.columns = pop_dat.columns.str.lower()


for v in ['intime','outtime','admittime']:
       pop_dat[v] = pd.to_datetime(pop_dat[v])

# do any additional filtering beyond cohort construction, if desired...

print('cohort loaded')
print("n")
print(len(final_ICU_IDs))


###########
########### all vital/lab values have been loaded; preprocessing of most extreme values
###########


VERBOSE=1
for v in all_ts_vars:
	dat = all_ts_dats[v]
	if 'GCS' not in v: 
		if VERBOSE:
			print(v)

		vals = np.array(dat['value'],"float")
		LOWER_Q = 1
		UPPER_Q = 99
		lq,uq = np.percentile(vals,[LOWER_Q,UPPER_Q])
		if VERBOSE:
			print(lq,uq)
		dat.loc[vals<=lq,'value'] = lq
		dat.loc[vals>=uq,'value'] = uq

		#some manual clips for certain vars
		if v=='fio2':
			dat.loc[vals<=1,'value'] *= 100 #units off
			dat.loc[vals<21,'value'] = 21 #implausible (and rarely seen) to be lower than room air
		if v=='map':
			dat.loc[vals<=40,'value'] = 40 #implausible
		if v=='sbp':
			dat.loc[vals<=60,'value'] = 60 #implausible
		if v=='dbp':
			dat.loc[vals<=30,'value'] = 30 #implausible		
		if v=='gfr':
			dat.loc[vals>=100,'value'] = 100 #GFR doesn't really mean anything when that high		



###################
################### load in & process all actions. NOTE that times vary across trajectories now!
###################

def total_pressor_normed_amount(start_t,end_t,v_starts,v_ends,v_rates):
	# get total amount of pressor in specified period [t_start,t_end]
	# normalized by size of period
	# TODO: vectorize & speedup...?

	delta_t = end_t - start_t #+ 1e-50

	#integrate total pressors given this period, catching all edge cases
	pressor_amt = 0 #+ 1e-50
	for v_s, v_e, r in zip(v_starts,v_ends,v_rates):
		if v_s >= start_t and v_e <= end_t:
			pressor_amt += r*60*(v_e-v_s)		
		if v_s < start_t and v_e > end_t:
			pressor_amt += r*60*delta_t
		if v_s < start_t and v_e > start_t and v_e <= end_t:
			pressor_amt += r*60*(v_e-start_t)
		if v_s >= start_t and v_s < end_t and v_e > end_t:
			pressor_amt += r*60*(end_t-v_s)

	pressor_normed_amt = 1/delta_t*pressor_amt
	return pressor_normed_amt

#load in data; actions are already sorted by icu_id & start_time
#anonymous: double checked if sorted. yes
fluids_dat = pd.read_csv(PATH_TO_QUERY_DATA+'allfluids_and_bloodproducts_mv.csv')
fluids_dat['STARTTIME'] = pd.to_datetime(fluids_dat['STARTTIME'])
fluids_dat['ENDTIME'] = pd.to_datetime(fluids_dat['ENDTIME'])



# note that these starts and ends are defined according to final_ICU_IDs
# So if final_ICU_IDs is cut down later because of map not hypo
# or too few maps, then it will lead to a mismatch? yes bc code below
# uses indices

# will get the index of the first occurence of the id 
map_starts = map_dat['icustay_id'].searchsorted(final_ICU_IDs,'left')
map_ends = map_dat['icustay_id'].searchsorted(final_ICU_IDs,'right')

#filter & threshold
fluids_dat['AMOUNT'] += 1 #rounding edge cases
fluids_dat = fluids_dat.loc[fluids_dat['AMOUNT'] >= 250, :] #cut tiny fluids
#TODO: why is pandas throwing warning?
fluids_dat['AMOUNT'][fluids_dat['AMOUNT']>=2000] = 2000


vaso_dat = pd.read_csv(path_to_query_data+'vasopressors_mv.csv')
vaso_dat['STARTTIME'] = pd.to_datetime(vaso_dat['STARTTIME'])
vaso_dat['ENDTIME'] = pd.to_datetime(vaso_dat['ENDTIME'])

#threshold; don't filter bc *any* pressor on is a big deal, unlike fluids
vaso_dat['RATE_NORMED_NOREPI'][vaso_dat['RATE_NORMED_NOREPI']>=2.5] = 2.5

#cache for speedup
fluid_starts = fluids_dat['ICUSTAY_ID'].searchsorted(final_ICU_IDs,'left')
fluid_ends = fluids_dat['ICUSTAY_ID'].searchsorted(final_ICU_IDs,'right')
vaso_starts = vaso_dat['ICUSTAY_ID'].searchsorted(final_ICU_IDs,'left')
vaso_ends = vaso_dat['ICUSTAY_ID'].searchsorted(final_ICU_IDs,'right')

def get_raw_data_summary(id_check):
	rawmap=map_dat.loc[map_dat.icustay_id==id_check,]
	rawmap.iloc[1:20]
	rm=rawmap.set_index('charttime')
	rawvaso=vaso_dat.loc[vaso_dat.ICUSTAY_ID==id_check,].copy()
	rawvaso.set_index('STARTTIME',inplace=True)
	rawvasoend=vaso_dat.loc[vaso_dat.ICUSTAY_ID==id_check,].copy()
	rawvasoend.ENDTIME  = rawvasoend.ENDTIME - timedelta(milliseconds=300)
	rawvasoend.set_index('ENDTIME',inplace=True)

	#rv = rawvaso
	rawcoh=cohort_dat.loc[cohort_dat.ICUSTAY_ID==id_check,]
	rmm=rm.join(rawvaso,how='outer')#.loc[:,['valuenum','RATE_NORMED_NOREPI','ENDTIME']].iloc[0:20]
	#rm.join(rm,rawvasoend,how='outer').loc[:,['valuenum','RATE_NORMED_NOREPI']].iloc[0:20]
	rmse = rmm.join(rawvasoend,how='outer',lsuffix='l')
	subsrmse=rmse.loc[:,['valuenum','RATE_NORMED_NOREPI','RATE_NORMED_NOREPIl']].iloc[0:50]
	subsrmse.RATE_NORMED_NOREPI[np.isnan(subsrmse.RATE_NORMED_NOREPI)] = 0  
	subsrmse.RATE_NORMED_NOREPIl[np.isnan(subsrmse.RATE_NORMED_NOREPIl)] = 0  
	subsrmse['vasorate'] = subsrmse.RATE_NORMED_NOREPI+subsrmse.RATE_NORMED_NOREPIl
	subsrmse.columns = ['map','end_vaso','start_vaso','vasorate']
	subsrmse.to_csv(path_to_model_data+'exampleRawData201098.csv')
	pickle.dump(subsrmse,open(path_to_model_data+'exampleRawData.p','wb'))
	#if np.sum(subsrmse.vasorate)>0: pdb.set_trace()
	return subsrmse


get_raw_data_summary(id_check)


#first pass we go through and grab all the raw action amounts & times, then discretize
all_rawaction_data = {}


SICK_TIME_BUFFER = 1.01  

LONG_GAP_TIME = 4.01

loop_t = time()
for ID_ind,ID in enumerate(final_ICU_IDs):
	if ID_ind % 100 == 99:
		print("processing %d/%d, took %.2f" %(ID_ind+1,len(final_ICU_IDs),time()-loop_t))
		loop_t = time()
	# get ICU start time, out time, & LOS 
	# anonymous:don't like index by ind. why not ID. note need
	# pop dat to then match final_ICU_IDs
	this_pop_dat = pop_dat.iloc[ID_ind] 

	# rather than take entry into ICU, start at first hypotensive episode

	s = map_starts[ID_ind]
	e = map_ends[ID_ind]
	this_maps = map_dat[s:e]

	if hypo_start==1:
		start_time,_ = get_first_hypo(this_maps)
	else:
		start_time = this_pop_dat['intime']

	outtime = this_pop_dat['outtime']


	total_time = (outtime-start_time).total_seconds()/60/60

	#either discharge from ICU or capped length	
	# NOTE take off some time at very end to give us a buffer,
	# and allow us to assess the effect of the last action taken...
	# we will artifically force last decision time to be here, at very latest
	end_time = min(total_time,LOS_CAP)

	#get MAP values & times to help inform grid...
	map_times = np.array((this_maps['charttime'] - start_time).astype('timedelta64[m]').astype(float))/60
	map_vals = np.array(this_maps['value'])
	

	#filter to before end...
	map_inds = np.logical_and(map_times <= end_time, map_times >= 0)
	map_times = map_times[map_inds]
	map_vals = map_vals[map_inds]

	#UPDATE our terminating time so that last MAP is where we cut things off...
	# rather than end of stay (anonymous thinks)
	end_time = np.max(map_times)

	##### get treatment info for this patient

	#fluids
	s = fluid_starts[ID_ind]
	e = fluid_ends[ID_ind]
	this_fdat = fluids_dat[s:e]

	f_starts = np.array((this_fdat['STARTTIME'] - start_time).astype('timedelta64[m]').astype(float))/60
	f_amounts = np.array(this_fdat['AMOUNT'])

	#filter irrelevant fluids
	f_ind = f_starts < end_time
	f_starts = f_starts[f_ind]
	f_amounts = f_amounts[f_ind]

	#pressors

	s = vaso_starts[ID_ind]
	e = vaso_ends[ID_ind]
	this_vdat = vaso_dat[s:e]
	
	v_starts = np.array((this_vdat['STARTTIME'] - start_time).astype('timedelta64[m]').astype(float))/60
	v_ends = np.array((this_vdat['ENDTIME'] - start_time).astype('timedelta64[m]').astype(float))/60
	v_rates = np.array(this_vdat['RATE_NORMED_NOREPI'])


	#filter irrelevant pressors
	v_ind = v_starts < end_time
	v_starts = v_starts[v_ind]
	v_ends = v_ends[v_ind]
	v_rates = v_rates[v_ind]
	v_ends[v_ends > end_time] = end_time #also force pressors to end, at latest, at our artifical end

	### Step 1: we need to get decision times for building this trajectory,
	### 	start by getting all times when treatment decisions were made:
	###			- pressor started
	###			- pressor ended
	###			- fluid started (end irrelevant bc short)

	all_tx_times = np.unique(np.concatenate([f_starts,v_starts,v_ends])) 
	all_tx_times = all_tx_times[all_tx_times < end_time] #only allow tx actions before end buffer

	
	
	# I guess zero start time because we define start time as 0
	# hopung that negative times get converted to vaso before
	final_action_times = np.arange(2/300, end_time, TIME_GRID_INT_HRS).tolist()#range(0,int(end_time))
	#print(final_action_times)
	if len(final_action_times)==0: pdb.set_trace()
	final_action_times = [1/300] + final_action_times + [final_action_times[-1]+0.01]


	#if ID == 201098:	pdb.set_trace()
	### ok, now collect up the action amounts...

	#tricky so use helper func
	
	pressor_normed_amts = []
	for (t0,t1) in zip(final_action_times[:-1],final_action_times[1:]):
		amt = total_pressor_normed_amount(t0,t1,v_starts,v_ends,v_rates)
		pressor_normed_amts.append(amt)
	pressor_normed_amts.append(np.nan) #none at end; not tracking anything after terminal end_time

	#more straightforward; not integrating, just treat as point mass with no timing
	fluid_amts = []
	for (t0,t1) in zip(final_action_times[:-1],final_action_times[1:]):
		fluid_inds = np.logical_and(f_starts >= t0, f_starts < t1)
		fluid_amts.append(np.sum(f_amounts[fluid_inds]))
	fluid_amts.append(np.nan)

	############ 
	actions_dat = pd.DataFrame()
	actions_dat['Times'] = final_action_times
	actions_dat['Vasopressor_normed_amt'] = np.array(pressor_normed_amts)
	actions_dat['Total_fluid_bolus_amt'] = np.array(fluid_amts)

	### SAVE
	all_rawaction_data[ID] = actions_dat

pickle.dump(all_rawaction_data,open(path_to_model_data+'all_raw_actions_times.p','wb'))
all_rawaction_data = pickle.load(open(path_to_model_data+'all_raw_actions_times.p','rb'))

# anonymous added this. else get keyerror, since some icuids are not in all_rawaction_data I think.
# this should be ok, since we subset pop data too. pop data used later.


final_ICU_IDs = np.sort(np.array(list(all_rawaction_data.keys())))
pop_dat = pop_dat.loc[np.in1d(pop_dat['icustay_id'],final_ICU_IDs),:]
assert final_ICU_IDs.shape[0] == pop_dat.shape[0]
print("new n after map/action filter")
print(len(final_ICU_IDs))


# ok, so only about 0.5% of time is there an action btw 0 & first flagged action time we have...
# ignore and these won't be modeled explicitly w RL but will at least inform states 
# when we aggregate some features based on past tx amounts given...
# TODO: what the , how/why is this edge case happening...

########## ok, now let's discretize actions into bins...

fluid_cutoffs = [0,500,1000,1e8] #4 bins
vaso_cutoffs = [0,8.1,21.58,1e8] #4 bins based on 33.3/66.7% quantiles of nonzeros

EPS = 1e-8
FLUIDS_BINS = np.array(fluid_cutoffs)+EPS; FLUIDS_BINS[0] = 0
VASO_BINS = np.array(vaso_cutoffs)+EPS; VASO_BINS[0] = 0

NUM_VASO_BINS = len(VASO_BINS) #4
NUM_FLUID_BINS = len(FLUIDS_BINS) #4
#in overall numbering, we'll have 0 = no action; 1 = dose 1 of fluids, 0 vaso; 
# 2 = dose 2 of fluids, 0 vaso; ... ; F = max fluids, 0 vaso; F+1 = 0 fluids, 1 vaso, ... 
NUM_ACTIONS = NUM_VASO_BINS*NUM_FLUID_BINS

def vaso_fluid_to_bin(fluids,vasos):
	f_ids = np.searchsorted(FLUIDS_BINS,fluids)
	v_ids = np.searchsorted(VASO_BINS,vasos)

	overall_ids = f_ids + NUM_FLUID_BINS*v_ids

	return f_ids,v_ids,overall_ids


##### ok group raw action amounts...almost done...

all_actions_disc_data = {}

loop_t = time()

for ID_ind,ID in enumerate(final_ICU_IDs):
	if ID_ind % 100 == 99:
		print("processing %d/%d, took %.2f" %(ID_ind+1,len(final_ICU_IDs),time()-loop_t))
		loop_t = time()

	dat = all_rawaction_data[ID]
        
	fluids = np.array(dat['Total_fluid_bolus_amt'])
	vasos = np.array(dat['Vasopressor_normed_amt'])

	f_ids,v_ids,overall_ids = vaso_fluid_to_bin(fluids,vasos)
	f_ids[-1] = v_ids[-1] = overall_ids[-1] = -999 #no last action

	dat['OVERALL_ACTION_ID'] = overall_ids
	dat['FLUID_ID'] = f_ids
	dat['VASO_ID'] = v_ids

	all_actions_disc_data[ID] = dat

pickle.dump(all_actions_disc_data,open(path_to_model_data+'all_disc-v4f4_actions_times.p','wb'))
# all_actions_disc_data = pickle.load(open(path_to_model_data+'all_disc-v4f4_actions_times.p','rb'))

#get overall action cts
all_actions_disc = []
for ID in final_ICU_IDs:
	acts = np.array(all_actions_disc_data[ID]['OVERALL_ACTION_ID'])[1:-1]
	all_actions_disc.append(acts)
all_actions_disc = np.concatenate(all_actions_disc)

action_cts = np.unique(all_actions_disc,return_counts=True)
for a,c in zip(action_cts[0],action_cts[1]):
	print(a,c,c/len(all_actions_disc)*100)



##### go back and add in extra vars: total fluids & pressors so far, and in past 8 hours

### TODO: THIS IS BROKEN & NOT REALLY 8 HOURS BC OF TIME GAPS...ITS LAST 8 TIME PTS...
## anonymous: note this.  I am not currently using this though

def rolling_sum(vec,k=8):
	res = np.cumsum(np.array(vec,"float"))
	res[k:] -= res[:-k]
	return res

loop_t = time()
for ID_ind,ID in enumerate(final_ICU_IDs):
	if ID_ind % 100==99:
		print("processing %d/%d, took %.2f" %(ID_ind+1,len(final_ICU_IDs),time()-loop_t))
		loop_t = time()

	actions_dat = all_actions_disc_data[ID]

	### off by 1 here; the prev/last vars should be 0 initially, and then drop the last one...
	z = np.zeros(actions_dat.shape[0])

	actions_dat['total_all_prev_vasos'] = z
	actions_dat['total_all_prev_vasos'][1:] = np.cumsum(np.array(actions_dat['Vasopressor_normed_amt'],"float"))[:-1]
	actions_dat['total_all_prev_fluids'] = z
	actions_dat['total_all_prev_fluids'][1:] = np.cumsum(np.array(actions_dat['Total_fluid_bolus_amt'],"float"))[:-1]
	actions_dat['total_last_8hrs_vasos'] = z
	actions_dat['total_last_8hrs_vasos'][1:] = rolling_sum(actions_dat['Vasopressor_normed_amt'],8)[:-1]
	actions_dat['total_last_8hrs_fluids'] = z
	actions_dat['total_last_8hrs_fluids'][1:] = rolling_sum(actions_dat['Total_fluid_bolus_amt'],8)[:-1]

	all_actions_disc_data[ID] = actions_dat

pickle.dump(all_actions_disc_data,open(path_to_model_data+'all_disc-v4f4_actions_times_prevactions.p','wb'))
# all_actions_disc_data = pickle.load(open(path_to_model_data+'all_disc-v4f4_actions_times_prevactions.p','rb'))


########## yay lets do states already


all_state_data = {}

baseline_cov_names = ['age','is_F','surg_ICU','is_not_white',
	'is_emergency','is_urgent','hrs_from_admit_to_icu','died_in_icu','first_careunit','los','los_hosp']

#cache cutpoints in sorted dataframes in advance, much faster
starts_ts = {}; ends_ts = {}
for v in all_ts_vars:
	dat = all_ts_dats[v]
	starts_ts[v] = dat['icustay_id'].searchsorted(final_ICU_IDs,'left')
	ends_ts[v] = dat['icustay_id'].searchsorted(final_ICU_IDs,'right')

# keeps track of index where first low MAP is observed
first_low_maps = {}

# Joe's loop
loop_t = time()
for ID_ind,ID in enumerate(final_ICU_IDs):
	if ID_ind % 100==99:
		print("processing %d/%d, took %.2f" %(ID_ind+1,len(final_ICU_IDs),time()-loop_t))
		loop_t = time()

	###############
	#static vars first
	# note names in lower case, even tho upper in cohort data. just was convention here
	this_pop_dat = pop_dat.iloc[ID_ind]

	start_time = this_pop_dat['intime'] 
	# redefined start_time below to be first hypotensive state
	# however, use this start_time to get length of stay, eg
	outtime = this_pop_dat['outtime']
	total_time = (outtime-start_time).total_seconds()/60/60
	died_in_icu = this_pop_dat['died_in_icu']	
	first_careunit = this_pop_dat['first_careunit']
	
	#build out all static vars
	age = float(this_pop_dat['age'])
	if age > 90: age = 90 #BUG: tons of ages are 300!! will screw up standardization
	if np.isnan(age): age = base_medians['age']

	is_F = int(this_pop_dat['gender']=='F')
	if np.isnan(is_F): is_F = 0

	icu_type = this_pop_dat['first_careunit']
	surg_ICU = int(icu_type=='CSRU' or icu_type=='SICU' or icu_type=='TSICU')

	ethn = this_pop_dat['ethnicity']
	is_not_white = 0
	if 'WHITE' not in ethn: is_not_white = 1

	admit_type = this_pop_dat['admission_type']
	los = this_pop_dat['los']
	los_hosp  = this_pop_dat['los_hosp']
	is_emergency = 0
	is_urgent = 0
	
	# anonymous commented: 
       # Emergency: The patient requires immediate medical intervention as a result of severe, 
       # life threatening or potentially disabling conditions. Generally, the patient is admitted through the emergency room. 
	# Urgent: The patient requires immediate attention for the care and treatment of a physical or mental disorder.
	
	if admit_type=='EMERGENCY': is_emergency=1
	if admit_type=='URGENT': is_urgent=1
	
	

	hrs_from_admit_to_icu = (this_pop_dat['intime']-this_pop_dat['admittime']).total_seconds()/60/60
	if hrs_from_admit_to_icu<0: hrs_from_admit_to_icu = 0 #weird edge case

	baseline_covs = np.array([age,is_F,surg_ICU,
		is_not_white,is_emergency,is_urgent,hrs_from_admit_to_icu,died_in_icu,first_careunit,los,los_hosp])

	###############
	##### figure out time stamps

	this_act_dat = all_actions_disc_data[ID]

	#drop initial action time; just so we can get (rare) initial actions
	# anonymous, stopped dropping initial action time. we want bp there
	grid_times = np.array(this_act_dat['Times'])#[1:] 
	n_t = len(grid_times)
	n_dec = n_t - 1 #no action taken at very last time point; used for final transition & reward


	###############
	### state variables
	states_dat = pd.DataFrame()
	states_dat['Times'] = grid_times
	states_dat['normed_time'] = grid_times/LOS_CAP

	#add baseline data in first
	for var,val in zip(baseline_cov_names,baseline_covs):
		states_dat[var] = val


	# anonymous is redefining start time to be at onset of hypotension here, 
	# as done above iwth actions
	# get MAP values & times to help inform grid...

	#if ID == 201098: pdb.set_trace()
	s = map_starts[ID_ind]
	e = map_ends[ID_ind]
	this_maps = map_dat[s:e]

	if hypo_start==1:
		start_time,_ = get_first_hypo(this_maps)

	#####
	#now build out the time series
	#also build out indicator vector at each time, noting whether var was imputed or not
	for v in all_ts_vars:
		s = starts_ts[v][ID_ind]
		e = ends_ts[v][ID_ind]
		this_dat = all_ts_dats[v][s:e]
		this_t = np.array((this_dat['charttime'] - start_time).astype('timedelta64[m]').astype(float))/60
	
		if 'GCS' in v:
			this_vals = np.array(this_dat['valuenum'])
		else:
			this_vals = np.array(this_dat['value'])

		#impute anything initially missing with pop median, then
		#fill this in with observed values via LOCF
		imputed_vals = ts_medians[v]*np.ones(n_t)

		### Now get LOCF for cts labs/vitals
		tt = np.searchsorted(grid_times+1e-8,this_t)
		# gives index of where each time where this_t goes into grid times
		
		# I think ignore this for map, etc. those dealt w below.
		# if multple values
		for i in range(len(tt)):
			if i!=len(tt)-1:
				imputed_vals[tt[i]:tt[i+1]] = this_vals[i] 
			else:
				imputed_vals[tt[i]:] = this_vals[i] 
				
		#EXCEPT: for MAP/SBP/DBP, we fill in with the worst in the window
		#	this only applies to windows in which *more* than 1 MAP is measured.
		#	after these windows, the LOCF kicks in with the most recent value from them

		# I wonder if actually taking -1 below amounts to LOCF above..
		# here we are just rewriting?


		if v in ['map','sbp','dbp']:
			#if (ID == 201098 and v=='map'): pdb.set_trace()
		
			u_tt = np.unique(tt)
			starts = np.searchsorted(tt,u_tt,'left')
			ends = np.searchsorted(tt,u_tt,'right')
			for t,s,e in zip(u_tt,starts,ends):
				if t>=0 and t<n_t: 
					# anonymous changed from min to most recent
					# for our purposes better
					# since we base actions on this, so most
					# recent is better because may have had a lot of
					# history before that not relevant to next action
					#  maybe mean is stabler. don't like min that much
					#if (ID==201098 and v=='map'):pdb.set_trace()
					#imputed_vals[t] = np.min(this_vals[s:e])
					#imputed_vals[t] = np.mean(this_vals[s:e])
					imputed_vals[t] = this_vals[s:e][-1]

		#get indicators for at which times the variable was actually sampled
		inds_samples_vals = np.zeros(n_t+1) #edge case when values past endtime
		inds_samples_vals[tt] = 1.0
		inds_samples_vals = inds_samples_vals[:-1] #edge case when values past endtime

		states_dat[v] = imputed_vals
		states_dat[v+'_ind'] = inds_samples_vals

	#last, combine the 3 GCS vars & then toss the individuals
	states_dat['GCS'] = states_dat['GCS_eye']+states_dat['GCS_motor']+states_dat['GCS_verbal']
	states_dat['GCS_ind'] = states_dat['GCS_eye_ind']+states_dat['GCS_motor_ind']+states_dat['GCS_verbal_ind']
	states_dat['GCS_ind'] = np.minimum(1,states_dat['GCS_ind'])
	
	#also drop GFR indicator since same as creat
	# anonymous commented out gfr_ind, since it does not appear to be there anymore...
	#states_dat = states_dat.drop(['GCS_eye','GCS_motor','GCS_verbal','GCS_eye_ind',
	#	'GCS_motor_ind','GCS_verbal_ind','gfr_ind'],axis=1)
	states_dat = states_dat.drop(['GCS_eye','GCS_motor','GCS_verbal','GCS_eye_ind',
		'GCS_motor_ind','GCS_verbal_ind'],axis=1)

	this_states_maps = states_dat.map
        # anonymous: get the ix of the first low MAP
        # we want to start the trajectory here further down the code, so just record it here
	ltSICKTHRESH = this_states_maps[this_states_maps<=MAP_SICK_THRESH]
	states_dat['first_lowMAPtime'] = np.NAN
	states_dat['first_lowMAPix'] = np.NAN
	states_dat['time_FirstLowMApToDischOrCap'] = np.NAN
	# we are filling a dictionary, so can just take patients with MAP measurements
	if not ltSICKTHRESH.empty:
                # this is cumbersome...try cleaner way
		first_low_ix=ltSICKTHRESH.index[0]
		#first_low_loc=this_maps.index.get_loc(first_low_ix)
		first_low_maps[ID] = first_low_ix
		time_of_firstlowmap=states_dat.Times[first_low_ix]
		states_dat['first_lowMAPtime'] = time_of_firstlowmap
		states_dat['first_lowMAPix'] = first_low_ix
		states_dat['time_FirstLowMApToDischOrCap'] = los - time_of_firstlowmap
	
	all_state_data[ID] = states_dat

###

###	

# THIS SHOULD NO LONGER CHANGE ANYTHING. WE DO SELECTION BEFORE STATE PROCESSING
print("Redefining cohort after state")
final_ICU_IDs = np.sort(np.array(list(all_state_data.keys())))
###update pop_dat to account for filtering from stays with few MAPs...
pop_dat = pop_dat.loc[np.in1d(pop_dat['icustay_id'],final_ICU_IDs),:]
assert final_ICU_IDs.shape[0] == pop_dat.shape[0]

pickle.dump(all_state_data,open(path_to_model_data+'all_states.p','wb'))
# all_state_data = pickle.load(open(path_to_model_data+'all_states.p','rb'))


############
############ setup rewards for each ICU stay
############

def lin_reward_func(bps,cutoffs=[40,55,60,65],vals=[-1,-0.15,-0.05,0]): 
	return np.interp(bps,cutoffs,vals)
# xx = np.linspace(40,75,1000); plt.plot(xx,np.log(1+lin_reward_func(xx)+1e-8)); plt.show()
# xx = np.linspace(40,75,1000); plt.plot(xx,lin_reward_func(xx)); plt.show()

all_reward_data = {}

URINE_OUTPUT_THRESH = 30 #if urine is above this, we're not too worried about BP if it's above 55
MAP_UO_THRESH = 55 #as long as MAP is above this, give max reward as long as UO is ok

for ID in final_ICU_IDs:

	s_dat = all_state_data[ID]

	times = np.array(s_dat['Times'])#[:-1] #last state doesn't have reward, no action (just use state to compute reward)
	
	# whether to onlc calc reward on bpt+1. 
	#Matters if want original scheme where r(s_{t+1}) is reward.
	# if rbptp1=1, then calc r(s_t)
	rbptp1 = 0
	if rbptp1:
		bps = np.array(s_dat['map'])[1:] #first bp reading not used for reward; bp_t, take a_t, r_t = f(bp_t+1)
	# anonymous is now using r_t=f(bp_t)
	else:
		bps = np.array(s_dat['map']) 
	

	#if ID == 201098:pdb.set_trace()
	if delta_rew:
		# we normalize by initial value, or else will try to act when map large
		# because will give large difference
		# doon t need to adjust for baseline map here
		#rewards = np.array([np.NAN] +list((bps[1:]-bps[:-1])/bps[:-1]))#[np.NAN]+bps[1:]-bps[:-1]
		# I think it s better to do  R = s'
		rewards = np.array([np.NAN] +list(bps[1:]))#[np.NAN]+bps[1:]-bps[:-1]
	elif last_map_rew:
		# need to adjust for baseline map here
		rewards = bps
	else:
		rewards = lin_reward_func(bps)
	

	#extra mask to ensure UO measured, if not measured yet, treat as if bad.
	# don't want to give a free pass to slightly low MAPs that we got
	# before a UO was measured

	if urine_in_reward:
		uos = np.array(s_dat['urine'])
		uo_inds = np.array(s_dat['urine_ind'])
		if np.all(uo_inds==0):
			uos *= 0
		else:
			first_meas_uo = np.where(uo_inds==1)[0][0]
			uos[:first_meas_uo] = 0
		#ok, now that fixed initial UOs, cut to one ahead for rewards...
		if rbptp1:
			uos = uos[1:]
		else:
			uos = uos

		good_inds = np.logical_and(uos >= URINE_OUTPUT_THRESH, bps >= MAP_UO_THRESH)

		rewards[good_inds] = 0 #moderate hypotension but good UO = ok

	rewards_dat = pd.DataFrame()
	rewards_dat['Times'] = times
	rewards_dat['Rewards'] = rewards

	all_reward_data[ID] = rewards_dat
idoi = 51
all_actions_disc_data[final_ICU_IDs[idoi]].loc[:,["Times","Vasopressor_normed_amt"]]
all_reward_data[final_ICU_IDs[idoi]]
all_state_data[final_ICU_IDs[idoi]].loc[:,['Times','map']]	

pickle.dump(all_reward_data,open(path_to_model_data+'rewards.p','wb'))
all_reward_data = pickle.load(open(path_to_model_data+'rewards.p','rb'))


##### extend the state space...
##### build out additional indicator variables for labs: 8 hour and ever-measured inds
#####	(no vitals since super frequent as is)


### helper funcs to convert inds at hourly level to other levels
def convert_1inds_to_kinds(ind_vec,times,k=8):
	#NOTE: old func breaks here bc times not hourly...
	n = len(ind_vec)
	assert len(times)==n

	res = np.zeros(n,"float")
	res[0] = ind_vec[0]
	for i in range(1,n):
		this_t = times[i]
		rel_inds = np.logical_and(times >= this_t-k, times <= this_t)
		res[i] = np.any(ind_vec[rel_inds]==1)

	return res

def convert_1inds_to_everinds(ind_vec):
	res = np.cumsum(np.array(ind_vec))
	return np.array(res>=1,"float")

#add in 8hour and ever-measured indicators for these variables (mostly labs, not measured as often)
#extra_ind_vars = np.array(['bicarbonate', 'bun', 'creatinine', 'fio2',
#       'glucose', 'hct', 'lactate', 'magnesium','platelets', 'potassium', 'sodium', 
#       'wbc', 'alt','ast', 'bilirubin_total', 'co2', 'hgb','inr','pco2', 'po2', 'weight'])

# anonymous removed co2. 

extra_ind_vars = np.array(['bicarbonate', 'bun', 'creatinine', 'fio2',
       'glucose', 'hct', 'lactate', 'magnesium','platelets', 'potassium', 'sodium',
       'wbc', 'alt','ast', 'bilirubin_total', #'co2', 
	'hgb','inr','pco2', 'po2', 'weight'])

all_extended_states_dat = {}
loop_t = time()
for ID_ind,ID in enumerate(final_ICU_IDs):
	if ID_ind % 100==99:
		print("processing %d/%d, took %.2f" %(ID_ind+1,len(final_ICU_IDs),time()-loop_t))
		loop_t = time()

	states_dat = all_state_data[ID]

	this_times = np.array(states_dat['Times'])
	for v in extra_ind_vars:
		# SLOW!!
		# states_dat[v+'_8ind'] = convert_1inds_to_kinds(states_dat[v+'_ind'],this_times,8)
		states_dat[v+'_everind'] = convert_1inds_to_everinds(states_dat[v+'_ind'])
	acts_dat = all_actions_disc_data[ID]

	# this may be wrong. i just took off -1
	# I think it looks correct, adds up right. still unsure - not really using this too much
	#TODO check for off by 1's
	for act_v in range(1,NUM_VASO_BINS):
		#states_dat['last_vaso_'+str(act_v)] = np.array(np.array(acts_dat['VASO_ID'])[:-1]==act_v,"float")
		states_dat['last_vaso_'+str(act_v)] = np.array(np.array(acts_dat['VASO_ID'])==act_v,"float")

	for act_f in range(1,NUM_FLUID_BINS):
		#states_dat['last_fluid_'+str(act_f)] = np.array(np.array(acts_dat['FLUID_ID'])[:-1]==act_f,"float")
		states_dat['last_fluid_'+str(act_f)] = np.array(np.array(acts_dat['FLUID_ID'])==act_f,"float")

	# states_dat['total_all_prev_vasos'] = np.array(acts_dat['total_all_prev_vasos'])[1:]
	# states_dat['total_all_prev_fluids'] = np.array(acts_dat['total_all_prev_fluids'])[1:]
	# states_dat['total_last_8hrs_vasos'] = np.array(acts_dat['total_last_8hrs_vasos'])[1:]
	# states_dat['total_last_8hrs_fluids'] = np.array(acts_dat['total_last_8hrs_fluids'])[1:]

	states_dat['total_all_prev_vasos'] = np.array(acts_dat['total_all_prev_vasos'])#[1:]
	states_dat['total_all_prev_fluids'] = np.array(acts_dat['total_all_prev_fluids'])#[1:]
	states_dat['total_last_8hrs_vasos'] = np.array(acts_dat['total_last_8hrs_vasos'])#[1:]
	states_dat['total_last_8hrs_fluids'] = np.array(acts_dat['total_last_8hrs_fluids'])#[1:]

	all_extended_states_dat[ID] = states_dat

#list of all vars
state_vars = np.array(all_extended_states_dat[final_ICU_IDs[0]].columns)

pickle.dump(all_extended_states_dat,open(path_to_model_data+'all_states_extravars.p','wb'))
pickle.dump(first_low_maps,open(path_to_model_data+'all_first_low_map.p','wb'))
all_state_data = pickle.load(open(path_to_model_data+'all_states_extravars.p','rb'))


#############
############# write out as npz for safe loading
#############

all_states_np = []
all_actions_np = []
all_rewards_np = []
all_times_np = []
IDs = []
for ID in final_ICU_IDs:

	this_s = all_extended_states_dat[ID]
	this_a = all_actions_disc_data[ID]
	this_r = all_reward_data[ID]
	this_t = all_rawaction_data[ID].Times.values

	states_np = np.array(this_s)#[:,1:] # dim: T x n_S
	if rbptp1:
		act_np = np.array(this_a)[1:-1,3] #dim: T
		# np.array(this_a)[:-1,3]
	else:
		act_np = np.array(this_a)[1:-1,3] #dim: T
	if rbptp1:
		rew_np = np.array(this_r)[:,1] # dim: T
	else:
		rew_np = np.array(this_r)[1:,1] # dim: T
	actt_np = np.array(this_t)[1:-1]
	#if len(states_np)>1:
		# only if some states occured

	all_states_np.append(states_np)
	
	all_actions_np.append(act_np)
	all_rewards_np.append(rew_np)
	all_times_np.append(actt_np)
	IDs.append(ID)


np.savez(path_to_model_data+'states_actions_rewards_IDs.npz',
	all_states=all_states_np,all_actions=all_actions_np,
	all_action_times=all_times_np,
	all_rewards=all_rewards_np,all_IDs=IDs,
	state_var_names=state_vars)

print("Completed")

pdb.set_trace()
#sanity check

# take an id which had vasos
idoi = 51
id_check = final_ICU_IDs[idoi]
#id_check = final_ICU_IDs[100]


# now look at the "raw" data (before this script)
rawmap=map_dat.loc[map_dat.icustay_id==id_check,]
rawmap.iloc[1:20]
rm=rawmap.set_index('charttime')
rawvaso=vaso_dat.loc[vaso_dat.ICUSTAY_ID==id_check,]
rawvaso.set_index('STARTTIME',inplace=True)
rv = rawvaso
rawcoh=cohort_dat.loc[cohort_dat.ICUSTAY_ID==id_check,]
rm.join(rawvaso,how='outer').loc[:,['valuenum','RATE_NORMED_NOREPI']].iloc[0:20]
pdb.set_trace()
#pd.DataFrame({'actions':acs,'states':sts,'rewards':ree.Rewards})

#now compare demographics
#stt.head(1)
#rawcoh
# age, sex matches

# now map 
#stt.map #starts at .25
rawmap.iloc[7:15]
# starting at hyp and going forward, matches

# vaso
rawvaso.iloc[0:5]
#acc



# the times seem off.  let's look at what is put into the aciton, states 
[el for el in zip(state_vars[1:],all_states_np[idoi][0])]


# These are lined up now such that we have the state from time t - 1
# the action at time t, and the reward from the state at time t
# (re imputation, note that in a time window with multiple maps, we take the last one)
traj_act = all_actions_np[idoi]
traj_rew = all_rewards_np[idoi]
traj_sta = [el[19] for el in all_states_np[idoi]]

ID = 201098
this_s = all_extended_states_dat[ID]
this_a = all_actions_disc_data[ID]
this_r = all_reward_data[ID]
this_t = all_rawaction_data[ID].Times.values
np.array(this_a)[1:-1,3] 
np.array(this_r)[1:,1] 
lin_reward_func([58,48,53,67])
# so after seeing 58, I take 8 and get 48 and r(48)=-0.55
# after seeing 48, I take 12 and get 53 and r(53)=-0.26
# after seeing 53, I take 12 again and get 67 and r(67)=
# so it's 58, 8, -0.55
	  #48, 12, -0.26

#pd.DataFrame({"actions":traj_act,"rewards":traj_rew,"states":traj_sta})

all_actions_disc_data[201098]

all_actions_disc_data[final_ICU_IDs[idoi]].loc[:,["Times","Vasopressor_normed_amt"]]
all_reward_data[final_ICU_IDs[idoi]]
all_state_data[final_ICU_IDs[idoi]].loc[:,['Times','map']]	


pdb.set_trace()
# anonymous added
lts = [len(el) for el in all_rewards_np]

#### log, then standardize appropriate columns after filtering 

all_states = []
for ID in final_ICU_IDs:
	this_states = all_extended_states_dat[ID]
	all_states.append(np.array(this_states))
all_states = np.vstack(all_states)

all_state_vars = np.array(this_states.columns)

vars_no_processing = ['Times','normed_time','is_F','surg_ICU','is_not_white','is_emergency','is_urgent','last_reward']
inds_to_process = [] #normed_time (already scaled); anything with ind; last_vaso_X or last_fluid_X;
for i,v in enumerate(all_state_vars):
	if v in vars_no_processing or 'ind' in v or 'last_fluid' in v or 'last_vaso' in v:
		pass
	else:
		inds_to_process.append(i)
		#print(all_state_vars[i])
inds_to_process = np.array(inds_to_process)

# all_state_means = np.mean(all_states,0)
# all_state_sds = np.std(all_states,0)
# for v,m,s in zip(all_state_vars,all_state_means,all_state_sds):
# 	print(v,m,s)

all_logstate_means = np.mean(np.log(all_states+.1),0)
all_logstate_sds = np.std(np.log(all_states+.1),0)


