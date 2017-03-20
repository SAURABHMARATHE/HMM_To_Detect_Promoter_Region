# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from hmmlearn import hmm
from sklearn.preprocessing import LabelEncoder
from random import randint

# S = Start State, E = Exit State, F = First Nucleotide State, N = Non Codon State, C = Codon State
states = ["S", "P", "N", "E"]
n_states = len(states)

emissions = ["A", "T", "C", "G", "N"]
n_emissions = len(emissions)
# Convert character emission symbols to integers as 
# required by MultinomialHMM
discrete_emissions = LabelEncoder().fit_transform(emissions)
print "LabelEncoder().fit_transform(", emissions, ") = ", discrete_emissions

# Create discrete emission dictionary
emission_dict = dict(zip(emissions, discrete_emissions))
print "emission_dict = ", emission_dict


start_probability = np.array([1.0])

transition_probability = np.array([
  [0.0, 0.50, 0.50, 0.0],
  [0.0, 0.50, 0.40, 0.10],
  [0.0, 0.40, 0.50,0.10], 
  [0.0, 0.0, 0.0, 0.0],
  
])

emission_probability = np.array([
  [0.20, 0.20, 0.20, 0.20, 0.20],
  [0.15, 0.15, 0.35, 0.35, 0.0],
  [0.35, 0.35, 0.10, 0.10, 0.10],
  [0.25, 0.25, 0.25, 0.25, 0.0],
  
])

choice=raw_input("Please enter choice(1-W/O training, 2-With training):")

if(choice=='1'):
    model = hmm.MultinomialHMM(n_components = n_states, params = "ste", n_iter = 0)
elif(choice=='2'):
    model = hmm.MultinomialHMM(n_components = n_states, params = "ste", init_params = "ste", n_iter = 50)
else:
    print "Wrong choice, Run prog again"
        
#model = hmm.MultinomialHMM(n_components = n_states, params = "ste", \
#        n_iter = 0)
#model = hmm.MultinomialHMM(n_components = n_states, params = "ste", \
#        init_params = "ste", n_iter = 50)
model.startprob_ = start_probability
model.transmat_ = transition_probability
model.emissionprob_ = emission_probability

# Read in sample sequence from file
sample_seq1 = []
delimiters=['\n','\r',' ']
with open('hmm_dataset.txt', 'r') as f: 
    # Read in entire sequence in file and remove all carriage returns "\n"
    sample_seq1 = list(f.read().replace('\n', ''))
    f.closed

sample_seq1=[x for x in sample_seq1 if x not in delimiters]

#print "sample_seq before conversion = ", sample_seq

# Convert sample seq to integer numbers per emission_dict
sample_seq2 = [emission_dict[e] for e in sample_seq1]
#print "sample_seq after conversion = ", sample_seq

sample_seq = np.array([sample_seq2]).T

"""
Estimate model parameters.
An initialization step is performed before entering the EM algorithm. 
If you want to avoid this step, set the keyword argument init_params 
to the empty string ‘’. Likewise, if you would like just to do an initialization, 
call this method with n_iter=0.
"""
model = model.fit(sample_seq)

print ">>> Get parameters for the estimator."
model_parameters = model.get_params(deep = True)
print "Model Parameters = ", model_parameters
print "Changed probabilities:",model.emissionprob_
print ""

logprob, hidden_model = model.decode(sample_seq, algorithm="viterbi")
print "Sample Sequence (Characters:", sample_seq1
print "Length Of Sample Sequence = ", len(sample_seq1)
print "Sample Sequence (Integers): ", sample_seq2
#print "Hidden Model (i.e., from  Viterbi Algorithm):", ", ".join(map(lambda x: states[x], hidden_model))
print

state_sequence=[]
state_sequence=map(lambda x: states[x], hidden_model)
print "Hidden Model (i.e., from  Viterbi Algorithm):", state_sequence
print "-------------------------------------------------------"
print state_sequence
print "Indices where CG pairs encountered and so Promoter region detected is:\n"
i=1
count=0
while i <= (len(state_sequence)-2):
    if (state_sequence[i]=='P' and state_sequence[i+1]=='P' and i+1<len(state_sequence)):
        str=sample_seq1[i-1]
        index=i
        while state_sequence[i+1]=='P':
            str=str+sample_seq1[i]
            i=i+1
        str=str+sample_seq1[i]+sample_seq1[i+1]
        print "start index:",index,"->",str,"->",i
        count=count+1
        i=i+1
    else:
        i=i+1
print ""
print "Total no of CpG pairs encountered in Viterbi algo:",count

print ">>> Compute the log probability under the model."
model_score = model.score(sample_seq)
print "Model Score = ", model_score
print

print ">>> Compute the model posteriors (using Forward & Backward algorithm)."
print ">>> Just last posterior probs are printed by default. "
print ">>> To print all posterior probs, change _PRINT_ALL_ in source code to True."
_PRINT_ALL_ = False
model_posterior_probs = model.predict_proba(sample_seq)
if _PRINT_ALL_:
    print "Model Posterior Probss = ", model_posterior_probs
else:
    print "Last Model Posterior Prob = ", model_posterior_probs[len(model_posterior_probs)-1]
print "Length of Posteriors = ", len(model_posterior_probs)
print

print ">>> Find most likely state sequence corresponding to observed emissions using the Baum-Welch algorithm."
model_seq = model.predict(sample_seq)
print "Model Most Likely State Seq = ", model_seq

'''
print "CpG sequence as per Baum Welch state sequence is found at index:\n"
for i in range(1,len(model_seq)-1):
    if(model_seq[i]==1 and model_seq[i+1]==1 and model_seq[i-1]!=1):
        print i

print
'''
print "CpG sequence as per Baum Welch state sequence is found at index:\n"
i=1
count=0
while i <= (len(model_seq)-2):
    if (model_seq[i]==1 and model_seq[i+1]==1 and i+1<len(model_seq)):
        str=sample_seq1[i-1]
        index=i
        while model_seq[i+1]==1:
            str=str+sample_seq1[i]
            i=i+1
        str=str+sample_seq1[i]+sample_seq1[i+1]
        print "start index:",index,"->",str,"->",i
        count=count+1
        i=i+1
    else:
        i=i+1
print ""
print "Total no of CpG pairs encountered in Baum Welch algo:",count

print""
print "Model Score = ", model_score

