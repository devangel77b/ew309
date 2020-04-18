% test_prob_nhits.m
% Test that prob_nhits seems to be working as expected for simple cases
% D Evangelista (2020)
% evangeli at usna dot edu

%% simple test
prob_nhits(1,1,0.5) %== 0.5 
prob_nhits(0,1,0.5) %== 0.5
prob_nhits(1,1,0.1) %== 0.1
prob_nhits(0,1,0.1) %== 0.9
prob_nhits(1,1,0.9) %== 0.9
prob_nhits(0,1,0.9) %== 0.1 

%% another simple test
prob_nhits(1,2,0.5) %== (2/4) = 0.5
prob_nhits(0,2,0.5) %== (1/4) = 0.25
prob_nhits(2,2,0.5) %== (1/4) = 0.25
prob_nhits(1,2,0.1) %== 2*0.1*0.9 = 0.18
prob_nhits(0,2,0.1) %== 1*0.9*0.9 = 0.81
prob_nhits(2,2,0.1) %== 1*0.1*0.1 = 0.01
% miss miss
% hit miss
% miss hit
% hit hit ==> gives 2/4 have 1 hit = 0.5 

%% another simple test
prob_nhits(0,3,0.5) %== (1/8) = 0.125
prob_nhits(1,3,0.5) %== (3/8) = 0.375
prob_nhits(2,3,0.5) %== (3/8) = 0.375
prob_nhits(3,3,0.5) %== (1/8) = 0.125
% miss miss miss
% hit  miss miss
% miss hit  miss
% miss miss hit
% hit  hit  miss
% miss hit  hit 
% hit  miss hit 
% hit  hit  hit

%% check that things sum to 1
psum = 0;
for i=0:10
  psum = psum+prob_nhits(i,10,0.5);
end
psum % == 1.0

% check that things sum to 1
psum = 0;
for i=0:20
  psum = psum+prob_nhits(i,20,0.1);
end
psum % == 1.0

% check that things sum to 1
psum = 0;
for i=0:20
  psum = psum+prob_nhits(i,20,0.9);
end
psum %== 1.0