# DENGUE_CHINA_2019

Code for Haobo Ni, Xiaoyan Cai, Jiarong Ren, Tingting Dai, Jiayi Zhou, Jiumin Lin, Li Wang, Lingxi Wang, Yunchong Yao, Ting Xu, Lina Xiao, Qiyong Liu, Xiaobo Liu, Pi Guo, "Epidemiological characteristics and transmission dynamics of dengue fever in China: a nationwide modelling study", 2024.

Code to run inference for dengue fever spread in China at city level.

Special thanks to Prof. Sen Pei(http://www.columbia.edu/~sp3449/ and https://github.com/SenPei-CU ) for his help and improvements to the code!!

## Functions
1. infer.m: the main function to run inference
2. initialize.m: compute initial conditions for all subpopulations
3. initializepara_eakf.m: initialize the ensemble of parameters for EAKF
4. seeding.m: set up initial dengue infection for all subpopulations
5. model_eakf.cpp: run the transmission model
6. checkbound_para.m: resample values in ensemble members that are out of limits
7. checkbound.m: ensure that variables remain non-negative
8. lhsu.m: latin hypercube sampling

## Data and data structure
1. commutingdata.mat:
a. The commuting network structure is stored in two vectors: nl (neighbor list) and part (partition). The vector nl records all neighbors of cities connected in the commuting network, and the vector part specifies the index range for the neighbors of each city. To find the list of cities where commuters living in city i work, use nl(part(i):part(i+1)-1). For computational convenience, the first city in the list is city i itself.
b. The number of commuters is stored in the vector C. For instance, the number of commuters from city i to city j=nl(k) [where part(i)<=k<=part(i+1)-1] is C(k). Note the index k should be within the index range for city i, i.e., from part(i) to part(i+1)-1.
c. The average number of commuters between two cities is stored in the vector Cave. For instance, the average number of commuters between city i and city j=nl(k) [where part(i)<=k<=part(i+1)-1] is Cave(k). It is the average number of commuters in both directions.
2. population.mat: city population
3. dailyincidence.mat: the number of reported cases in each city on each day. 
4. parameter.mat: identified parameters (mosquito population birth rate, bite rate, the probability for transmission of the infection from mosquito to human, the probability for transmission of the infection from human to mosquito)
5. mosmun.mat: weighting of initial mosquito population in each city
6. weight.mat: initial parameter weights

## How to run inference
To speed up the inference, the transmission model is programmed in C++. The inference function calls the C++ function through interface between MATLAB and C++.

1. Compile the C++ function in MATLAB using “mex model_eakf.cpp”. A C++ compiler (e.g., Xcode) needs to be installed on the computer before running this command.
2. Run infer.m

## Outputs
1.inference.mat: model inference results for all parameters and variables
a. para_post: posterior parameters
b. Sm_post: posterior susceptible mosquito population in each city on each day
c. Im_post: posterior infected mosquito population in each city on each day
d. Sh_post: posterior susceptible human population in each city on each day
e. Ihr_post: posterior reported infected human population in each city on each day
f. Ihu_post: posterior unreported infected human population in each city on each day
g. dailyIm_post_rec: posterior estimate of the daily newly infected mosquito population in each city on each day
h. dailyIhr_post_rec: posterior estimate of the daily newly reported infected human population in each city on each day
i. dailyIhu_post_rec: posterior estimate of the daily newly unreported infected human population in each city on each day
j. dailyIm_prior_rec: prior estimate of the daily newly infected mosquito population in each city on each day
k. dailyIhr_prior_rec: prior estimate of the daily newly reported infected human population in each city on each day
l. dailyIhu_prior_rec: prior estimate of the daily newly unreported infected human population in each city on each day
m. dailyIhr_post_rec1: posterior estimates of newly reported infections per day by city as a result of population movements
n. dailyIhu_post_rec1: posterior estimates of newly unreported infections per day by city as a result of population movements
