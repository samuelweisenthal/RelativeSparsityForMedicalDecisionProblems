# This is a post processing file for the output of data_clean.py to put the python objects into form that
# is easily read by R.
# It takes the output from Futoma 2020, who wrote data_clean.py and puts it into a form
# that the script mimic.R can read for modeling.
# Also, we binarize action

import numpy as np
import pickle
import pandas as pd
import pdb
import os
import shutil

path = '/Users/anonymous/Box/MIMIC/mimic-iii-v1-4/hypotension-RL/model-data2'
os.chdir(path)
exp_name = ""
fluidAndVaso = 0

# creates multistage episodes

def make_multistage(d, name, maxT, at):
    # which states are floats? useful for putting states into np arrays, which only take floats/ints
    isfloatorint = [(isinstance(el, float) or isinstance(el, int)) for el in d['all_states'][0][0]]
    float_names = np.array([el for el in d['state_var_names']])[isfloatorint]
    pd.DataFrame({"float_names":float_names}).to_csv("float_names.csv")#,index=False)
    # start with hypotensive measurement
    states=[]

    actions = []
    rewards = []
    ids = []

    print("JUST TAKING FIRST ",maxT," DECISION POINTS")

    for (ss,aa,rr,id) in zip(d['all_states'],d['all_actions'],d['all_rewards'],d['all_IDs']):
        # only if at least 2 state measurements

        states.append(ss[0:][0:maxT])
     
        if fluidAndVaso:
            actions.append([1*(el>0) for el in aa[0:][0:maxT]])
        else:
            actions.append([1*(el>3) for el in aa[0:][0:maxT]])
        rewards.append(rr[0:][0:maxT])
        ids.append(id)

    pd.DataFrame(actions).to_csv((name+"multistage_actions.csv"))#,index=False)
    pd.DataFrame(rewards).to_csv((name+"multistage_rewards.csv"))#,index=False)
  
    #atkeys = [el for el in at.keys()]
    atkeys = d['all_IDs']

    #assert (atkeys == d['all_IDs'])
    Traj_times = at#[at[i] for i in atkeys]
    times = []
    vasoamts = []
    vaso_norm_amounts = []
    for atel in at:
        times.append(atel[0:maxT])
        #vasoamts.append(at[id].Vasopressor_normed_amt.values[0:maxT])

    pd.DataFrame(times).to_csv((name+"multistage_times.csv"))#,index=False)
    #pd.DataFrame(vasoamts).to_csv((name+"multistage_vasoamts.csv"))#,index=False)

   
    # trying to put trajectories into an object R can read
    # ended up instead just sending each time step to a different csv
    # send each time step to different csv (avoid reading numpy arrays into R)
    # todo: not just floats anymore


    lens=[len(el) for el in states]
    N=len(states)
    T = max(lens)
    K =  states[0][0].shape[0]
    states_floats = []

    for n in range(0,N):
        states_floats_n = []
        for t in range(0,len(states[n])):
            states_floats_n.append(states[n][t][isfloatorint])
        states_floats.append(states_floats_n)

    N=len(states_floats)
    T = max(lens)
    K =  states_floats[0][0].shape[0]
    ne = np.empty((N,T,K))
    ne[:] = np.NAN

    for n in range(0,N):
        for t in range(0,len(states_floats[n])):
            for k in range(0,K):
                ne[n,t,k] = states_floats[n][t][k]


    newpath = r'stages'+name 
    if os.path.exists(newpath): shutil.rmtree(newpath)
    os.makedirs(newpath)
    for t in range(0,maxT):
        tdf = pd.DataFrame(ne[:,t,:])
        tdf.columns = float_names
        time_name = "./"+newpath+"/"+str(t) + ".csv"
        tdf.to_csv(time_name)#,index=False)
    print("ms saved ok")
    return [ids,actions,rewards,ne]


# maximum number of decision points to take in the trajectories starting at hypotension onset
# should be number of states time steps - 1 , or number of action or reward time steps
maxT = 4


d=np.load('states_actions_rewards_IDs.npz',allow_pickle=True)
idoi = 51
# can use the id below to check
id_check = [el for el in d['all_IDs']][idoi]
mos=[el for el in d['all_states'][idoi]]
mos=[el[19] for el in mos]




moa=[el for el in d['all_actions'][idoi]]
mor=[el for el in d['all_rewards'][idoi]]



first_low_maps = pickle.load(open('all_first_low_map.p','rb'))
s=pickle.load(open('all_states.p','rb')) #used?
print([el for el in d.keys()])

#no_ind = [(('ind' not in el) and ('last_' not in el)) for el in d['state_var_names']]
# by data_clean.py, actions 0-3 are 0 vaso
# hence we need to check for each visit whether there is a number greater than 3
## why some zero. if no action.
print("n",len(d['all_actions']))

#this should check for emptiness, not take max
mxs=[]
for el in d['all_actions']:
	try:
		mxs.append(np.max(el))
	except: 
		mxs.append(np.NAN)

takes = [~np.isnan(el) for el in mxs]
did_sub=d['all_IDs'][takes]
dst_sub=d['all_states'][takes]
dac_sub=d['all_actions'][takes]
dre_sub=d['all_rewards'][takes]
dti_sub=d['all_action_times'][takes]

[el for el in d['all_IDs']][42]

d={'all_states':dst_sub,'all_actions':dac_sub,'all_rewards':dre_sub,
'all_IDs':did_sub,'all_action_times':dti_sub,'state_var_names':d['state_var_names']}
print("n",len(d['all_actions']))

idoi = 19
# can use the id below to check
id_check = [el for el in d['all_IDs']][idoi]
mos=[el for el in d['all_states'][idoi]]
mos=[el[19] for el in mos]
moa=[el for el in d['all_actions'][idoi]]
mor=[el for el in d['all_rewards'][idoi]]

# now 201098 is index 49, which would be index 50 in R 
#[el for el in d['all_IDs']][49]

# if greater than three, implies vaso given
# if greater than zero, either vaso or fluid given


if fluidAndVaso:
    gt3=[np.max(el)>0 for el in d['all_actions']]
else:
    gt3=[np.max(el)>3 for el in d['all_actions']]
np.mean(gt3) #.27. So about 1/4 of time, use vaso. check with Joe's summary
# Joe has more vaso. It looks like maybe 3/4 vaso..

#fs = []
#for (sts,id) in zip(d['all_states'],d['all_IDs']):
#	pdb.set_trace()
#	fs.append(sts[first_low_maps[id]])

all_states_list = [el for el in d['all_states']]

fs = []
for (id,el) in zip(d['all_IDs'],d['all_states']):
    try:
        fs.append(el[0])
    except:
        fs.append([np.NAN])

#fs=[el[0] for el in d['all_states']]
#fs_noind=[el[no_ind] for el in fs]
#varsnames_noind=d['state_var_names'][no_ind]
varsnames = d['state_var_names']

#mean_rewards = [np.mean(el) for el in d['all_rewards']]

mean_rewards = []
for el in d['all_rewards']:
    try:    
        mean_rewards.append(np.mean(el))
    except:
        print("in except")
        mean_rewards.append(np.NAN)

np.mean(mean_rewards)
np.max(mean_rewards)
np.min(mean_rewards)

fs_arr = np.array(fs)
sdf = pd.DataFrame(fs)
sdf.columns = varsnames
sdf.to_csv("first_states.csv")#,index=False)

ls = []
for el in d['all_states']:
    try:
        ls.append(el[-1])
    except:
        ls.append(np.NAN)
 #[el[-1] for el in d['all_states']]

#ls_noind=[el[no_ind] for el in ls]

ls_arr = np.array(ls)
sdf = pd.DataFrame(ls)
sdf.columns = varsnames
sdf.to_csv("last_states.csv")#,index=False)

# don't set index to false, because if so R will remove blanks
pd.DataFrame({'r':mean_rewards}).to_csv('mean_rewards.csv')#,index=False)
gt3=[int(el) for el in gt3]
pd.DataFrame({'a':gt3}).to_csv('vaso.csv')#,index=False)


at=d['all_action_times']#pickle.load(open("all_raw_actions_times.p","rb"))
atkeys = d['all_IDs']#[el for el in at.keys()]

ms=make_multistage(d, exp_name,maxT=maxT,at=at)

Ts=[np.max(el) for el in at]

pd.DataFrame({'Ts':Ts}).to_csv('Ts.csv')#,index=False)

pdb.set_trace()
tgt24 = {}
for id in atkeys: tgt24[id] = at[id].Times>24


#fluid_starts = fluids_dat['ICUSTAY_ID'].searchsorted(final_ICU_IDs,'left')

# aggregate first and last 12 hours
# just one way to make a 2 stage decision
# I wonder if more stable
cp=1
if cp ==1:
    changepoints = {}
    for id in atkeys: 
        changepoints[id] = tgt24[id].searchsorted('True','left')
    #tgt24[atkeys[2]].searchsorted('True','left')
    states=[el for el in d['all_states']]
    ids=[el for el in d['all_IDs']]
    rewards = [el for el in d['all_rewards']]
    actions = [el for el in d['all_actions']]

    ss = list()
    for (id,state) in zip(ids,states):
        if changepoints[id] < len(state):
            ss.append(state[changepoints[id]])
        else:
            ss.append(np.NAN)

    #sss = 

    #pd.DataFrame()

    rr = list()
    for (i,id) in enumerate(ids):
        tmask = tgt24[id][1:-1]
        rr.append(
            [np.mean(rewards[i][~tmask]),np.mean(rewards[i][tmask])]
            )

    aa = list()
    for (i,id) in enumerate(ids):
        tmask = tgt24[id][1:-1]
        aa.append(
            [np.sum(actions[i][~tmask]>3)>0,np.sum(actions[i][tmask]>3)>0]
            )    


    tgt24[ids[2]]





### end of multistage data generation
    
#before_sts =  pd.DataFrame(np.array(mean_before_hyps))

#before_sts.columns = float_names
#before_sts.to_csv("mean_before_hyp_floats.csv")#,index=False)

#np.save("hyp_s",ne)

#print("n (more than 2 states - problematic)",len(hyp_ids))

# first and last state (eg, if going to treat whole trajectory)
# as one decision

hyp_fs = [el[0] for el in hyp_states]
hyp_ls = [el[-1] for el in hyp_states]

hyp_mxs=[]
for el in hyp_actions:
    try:
        hyp_mxs.append(np.max(el))
    except: 
        #print("exception")
        hyp_mxs.append(np.NAN)

# whether ever given vaso

hyp_vaso=[el>3 for el in hyp_mxs]
pctv = np.mean(hyp_vaso)

# average reward

hyp_mean_rewards = [np.mean(el) for el in hyp_rewards]
hyp_first_rewards = []
for el in hyp_rewards:
  
    if len(el)>=1:
        hyp_first_rewards.append(el[0]) 
    else:

        hyp_first_rewards.append(np.NAN)
hyp_first_actions = []
for el in hyp_actions:
    if len(el)>=1:
        hyp_first_actions.append(el[0])
    else:
        hyp_first_actions.append(np.NAN)

# overall mean
print("prop w mean reward small",np.mean([el<0 for el in mean_rewards]))
print("prop given vaso",np.mean(gt3))

# hypo mean
hyp_first_vaso = [el>3 for el in hyp_first_actions]
pctrlz = np.mean([el<0 for el in hyp_mean_rewards])
print("hypo prop given vaso",pctv)
print("hypo prop r<0",pctrlz)

# hypo first rewards,actins
print("hypo prop w first reward small",np.mean([el<0 for el in hyp_first_rewards]))
print("hypo prop given vaso first",np.mean([el>3 for el in hyp_first_actions]))


# some rewards are nan
# ecause they only have one state after hypotension.
# so no second state to measure a reward
# we will therefore subset by those with at least 1 state after hypo

# actually, don't do this.  keep in form with multistage above

ls = [len(el) for el in hyp_states]
stlg1 = [el>1 for el in ls]


hyp_as=pd.DataFrame({"a":hyp_vaso})#[stlg1]
hyp_as.to_csv("hyp_a.csv")#,index=False)

hyp_first_as=pd.DataFrame({"a":hyp_first_vaso})#[stlg1]
hyp_first_as.to_csv("hyp_first_a.csv")#,index=False)


hyp_rs = pd.DataFrame({"r":hyp_mean_rewards})#[stlg1]
hyp_rs.to_csv("hyp_r.csv")#,index=False)

hyp_first_rs = pd.DataFrame({"r":hyp_first_rewards})#[stlg1]
hyp_first_rs.to_csv("hyp_first_r.csv")#,index=False)

#hyp_fs_noind=[el[no_ind] for el in hyp_fs]
fsdf=pd.DataFrame(hyp_fs)#[stlg1]
fsdf.columns = varsnames
fsdf.to_csv("hyp_fs.csv")#,index=False)
#hyp_ls_noind=[el[no_ind] for el in hyp_ls]
lsdf=pd.DataFrame(hyp_ls)#[stlg1]
lsdf.columns = varsnames
lsdf.to_csv("hyp_ls.csv")#,index=False)

pdb.set_trace()
print("NEED TO GO BACK AD DEAL WITH NULL ACTIONS")
pdb.set_trace()

