
import torch
import gpytorch
import pandas as pd
import numpy as np
from gpytorch.mlls import SumMarginalLogLikelihood
import matplotlib.pyplot as plt
# from matplotlib import cm
# import matplotlib.colors as colors
import scipy.stats as ss
# import seaborn as sns
import sys


# TODO make an interface to R package here

###################
### Data import ###
###################

# replace with the path to your data
# the input needs to be a table with columns:
#   - uniqueID : peptide-ID
#   - x : temperature
#   - y : fraction non-denaturated
#   - condition : treatment/control
# data = pd.read_csv("Example_peptides_for_fit.csv")
#data = pd.read_csv("~/Documents/Phd/Experiments/018/python/HT_ATP_all.csv")    # example data from the TPP NPARC package

def FitGPs_function(sm, path_list, trypticity):

    #common_path = "/Users/moni/Documents/PhD/Experiments/"
    path = sys.argv[1] + "0" + np.array2string(path_list["Experiment"][path_list["sm"] == sm].unique()[0]) + "/" + sm + "/python_import/"
    path_output = sys.argv[2] + "0" + np.array2string(path_list["Experiment"][path_list["sm"] == sm].unique()[0]) + "/" + sm + "/python_output/"

    #trypticity = "FT"
    filename = sm + "_" + trypticity + ".csv"
    data = pd.read_csv(path + filename)
    # data['y'] = np.log(data['y'] + 1)
    
    # Filer out peptides for which null and alternative model coincide or not exactly two different conditions exist (TODO to be extended to more than 2 conditions)
    peps_all = data['uniqueID'].unique()
    pep_with_2_cond = peps_all[[len(data[data['uniqueID'] == pep]['condition'].unique()) == 2 for pep in peps_all]]
    
    print("Filtered out", len(peps_all) - len(pep_with_2_cond), "peptide(s) where not exactly 2 different conditions were available")
    
    data = data[[rowpep in pep_with_2_cond for rowpep in data['uniqueID']]]
    peptides2test = data['uniqueID'].unique()
    conds = np.insert(data['condition'].unique(), 0, "full model")
    
    ########################
    ### Model definition ###
    ########################
    
    # define a GP model using exact inference with constant mean function and squared exponential (RBF) kernel
    class ExactGPModel(gpytorch.models.ExactGP):
        def __init__(self, train_x, train_y, likelihood):
            super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
            self.mean_module = gpytorch.means.ConstantMean() # posterior mean need not be zero, prior has high influence only for low noise values - e.g. actual noise as fixed noise
            self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel())
    
        def forward(self, x):
            mean_x = self.mean_module(x)
            covar_x = self.covar_module(x)
            return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)
    
    # initialize likelihood and model for each peptide/protein over all conditions and per condition
    # order of list is 1. peptide, 2. condition (full, condition_1, condition_2)
    # Model is the GP model defined above with inputs x_train (temperatures) and y_train (intensities), the likelihood is a Gaussian:
    # GaussianLikelihood() uses homoscedasitc residuals, FixedNoiseGaussianLikelihood() prespecified possibly different residual variances at each input point
    model_list = []
    l_list = []
    for pep in peptides2test:
        # filter data to peptide
        df = data[data['uniqueID'] == pep]
    
        # define joint model:
        temp = torch.as_tensor(np.asarray(df['x'])).double()
        intens = torch.as_tensor(np.asarray(df['y'])).double()
        lik = gpytorch.likelihoods.GaussianLikelihood()
        # noise_het = torch.tensor([intens[temp == tt].var() for tt in temp])
        # lik = gpytorch.likelihoods.FixedNoiseGaussianLikelihood(noise = noise_het, learn_additional_noise = False) # whitout additional noise it learns a very high scale for the GP kernel to compensate (due to the constant meant) - noisy interpolations
        l_list.append(lik)
        model_list.append(ExactGPModel(temp, intens, lik))
    
        # define models specific for each treatment:
        for cond in df['condition'].unique(): # TODO pass condH0 and condH1 (from the R package) here for general usage
            dfcond = df[df['condition'] == cond]
            temp = torch.as_tensor(np.asarray(dfcond['x'])).double()
            intens = torch.as_tensor(np.asarray(dfcond['y'])).double()
            lik = gpytorch.likelihoods.GaussianLikelihood()
            # noise_het = torch.tensor([intens[temp == tt].var() for tt in temp])
            # lik = gpytorch.likelihoods.FixedNoiseGaussianLikelihood(noise =  noise_het, learn_additional_noise = True)
            l_list.append(lik)
            model_list.append(ExactGPModel(temp, intens, lik))
    
    
    # dimensions
    n_models = len(model_list)
    n_cond = len(conds)
    n_pep = len(peptides2test)
    assert n_models == n_cond * n_pep, "Error in building the model list."
    
    ######################
    ### Model training ###
    ######################
    
    # Combine different models in an IndependentModelList for simultaneous fitting
    model = gpytorch.models.IndependentModelList(*model_list)
    likelihood = gpytorch.likelihoods.LikelihoodList(*l_list)
    
    # Find model hyperparameters
    model.train()
    likelihood.train()
    
    # Optimize the models using the Adam optimizer
    optimizer = torch.optim.Adam([
        {'params': model.parameters()},  # Includes all submodel and all likelihood parameters
    ], lr=0.1)  # learning rate
    
    # Train
    sum_mll = SumMarginalLogLikelihood(likelihood, model) # objective function given by sum marginal log-likelihood over all models
    training_iterations = 1000
   # counter = 0
    loss_hist = [np.inf] * training_iterations
    for i in range(training_iterations):
        optimizer.zero_grad()
        output = model(*model.train_inputs)
        loss = -sum_mll(output, model.train_targets)
        loss.backward()
        print('Iter %d/%d - Loss: %.6f' % (i + 1, training_iterations, loss.item()))
        optimizer.step()
        #loss = -sum_mll(output, model.train_targets)
        #print('Iter %d/%d - Loss: %.3f' % (i + 1, training_iterations, loss.item()))
        loss_hist[i] = loss.item()
        if i > 1 and abs(loss_hist[i] - loss_hist[i-1]) <= 0.0000011:
            break
        
        
    plt.plot(range(training_iterations), loss_hist)
    #plt.show()
    plt.savefig(path_output + 'Figure1_' + trypticity + '.png')
    
    # Set into eval mode
    model.eval()
    likelihood.eval()
    
    # Show fitted parameters for one example model (note that these are transformed for the actual values):
    for param_name, param in model.models[0].named_parameters():
        print(f'Parameter name: {param_name:42} value = {param.item()}')
        # model.models[0].state_dict()  # state dictionary contains all traininable parameters
    
    
    ###########################################
    ### result data frame based on LRT test ###
    ###########################################
    
    # Calculate marginal likelihoods for each model
    mll_eval = SumMarginalLogLikelihood(likelihood, model)
    outputs = model(*model.train_inputs) # returns per model a MultivariateNormal containing the posterior mean and covariance
    targets = model.train_targets
    mlls_list =[mll(output, target).item() for mll, output, target in zip(mll_eval.mlls, outputs, targets)]
    
    # Calculate LRT test statistics
    df = 1 # degree of freedom for the the chi2-distribution in a likelihood-ratio test (LRT) - 1 is used e.g. spatialDE TODO find a suitable df
    MLL_all = pd.DataFrame({'peptide' : np.repeat(peptides2test, n_cond), 'condition' : np.tile(conds, n_pep), 'MLL' : mlls_list})
    MLL_peptide =  MLL_all.pivot(index = 'peptide', columns = 'condition', values = 'MLL')
    MLL_peptide['LR'] = -2 * (2 * MLL_peptide['full model'] - MLL_peptide[conds].sum(axis = 1)) # likelihood ratio test statistics
    MLL_peptide['pval'] = 1 - ss.chi2.cdf(MLL_peptide['LR'], df) # p-value from a LRT test
    MLL_peptide['BF'] = 1 / np.exp(- MLL_peptide['LR'] / 2) # Bayes - Factor
    
    MLL_peptide.to_csv(path_output + "MLL_" + sm + "_" + trypticity + ".csv", sep=',', index=True, header=True)
    
    ##########################################
    ###  result data frame based on F-test ###
    ##########################################
    
    # Compare model fits based on residual sum of squares
    x_train_list = [submodel.train_inputs for submodel in model.models]
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        # Predictions for all outcomes as a list
        predictions = likelihood(*model(*x_train_list))
    RSS =[]
    for submodel, prediction in zip(model.models, predictions):
        tr_y = submodel.train_targets.detach().numpy()
        RSS.append(((prediction.mean - tr_y) ** 2).sum().item())
    RSS_all = pd.DataFrame({'peptide' : np.repeat(peptides2test, n_cond), 'condition' : np.tile(conds, n_pep), 'RSS' : RSS})
    RSS_pep = RSS_all.pivot(index = 'peptide', columns = 'condition', values = 'RSS')
    RSS_pep['RSS1'] = RSS_pep[conds].sum(axis = 1) - RSS_pep['full model']
    RSS_pep['RSS0'] = RSS_pep['full model']
    
    # histogram of residuals - Gaussian ok
    plt.hist((torch.cat([p.train_targets for p in model.models]) - torch.cat([p.mean for p in likelihood(*model(*x_train_list))])).detach().numpy())
    #plt.show()
    plt.savefig(path_output + 'Figure2_' + trypticity + '.png')

    
    # calculate p-value based on RSS (F-test)
    # degree of freedom can be determined as in NPARC by exporting the Fstat to the PARC package and use it instead of the spline's one
    # df-values below may not be meaningful
    df1 = 4 # TODO determine df - empirical approach
    df2 = 8 # TODO determine df - empirical approach
    RSS_pep['deltaRSS'] = (RSS_pep['RSS0'] -  RSS_pep['RSS1'])
    RSS_pep['Fstat'] = RSS_pep['deltaRSS'] / RSS_pep['RSS1']
    RSS_pep['pval'] =  1 - ss.f.cdf(RSS_pep['Fstat'], df1, df2)
    print(RSS_pep)
    
    # compare p-values from the two tests
    plt.figure()
    plt.scatter(MLL_peptide['pval'], RSS_pep['pval'])
    plt.ylabel('p-val (F-test)')
    plt.xlabel('pval (LRT)')
    #plt.show()
    plt.savefig(path_output + 'Figure3_' + trypticity + '.png')

    #################
    ### Plot fits ###
    #################
    
    # Make predictions (on selected test points - here a grid from minimum to maximum temperature)
    min_x = data['x'].min()
    max_x = data['x'].max()
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        # specify test point locations
        test_x = torch.linspace(min_x, max_x, 101).double()
        test_x_list = [test_x] * n_models # use the same x for all peptides and conditions
        # Predictions for all outcomes as a list
        flist = model(*test_x_list) # posterior process in each model
        fmean = torch.stack([f.mean for f in flist]) # should be the same as predictions_mean
        fconf_upper = torch.stack([f.confidence_region()[1] for f in flist])
        fconf_lower = torch.stack([f.confidence_region()[0] for f in flist])
        predictions = likelihood(*model(*test_x_list))
        predictions_mean = torch.stack([pred.mean for pred in predictions]) # array with fitted mean per peptide and condition
        predictions_conf_lower = torch.stack([pred.confidence_region()[0] for pred in predictions]) # array with lower confidence region per peptide and condition
        predictions_conf_upper = torch.stack([pred.confidence_region()[1] for pred in predictions]) # array with upper confidence region per peptide and condition
    
    # Plot samples for the first peptide and full model (i = 0)
    # Drawbacks: can sample negative values - might be useful to log-transform data for sampling purposes?
    i = 0
    plt.figure()
    plt.plot(test_x, fmean[i]) # mean
    plt.fill_between(test_x, fconf_upper[i], fconf_lower[i], alpha = 0.1) # mean
    plt.fill_between(test_x, predictions_conf_upper[i], predictions_conf_lower[i], alpha = 0.1) # mean
    for seed in range(5):
        plt.scatter(test_x, flist[i].sample(), s= 1) # samples from posterior GP
    for seed in range(5):
        plt.scatter(test_x, predictions[i].sample(), s=1)  # samples from marginal likelihood GP + sample noise
    #plt.show()
    plt.savefig(path_output + 'Figure4_' + trypticity + '.png')

    
    # data frame with fits (peptide
    df = pd.DataFrame({'peptide' : np.repeat(peptides2test, n_cond*len(test_x)),
                        'condition' : np.tile(np.repeat(conds, len(test_x)), n_pep),
                        't' : np.tile(test_x, n_pep * n_cond)})
    df['y'] = torch.flatten(predictions_mean)
    df['conf_lower'] = torch.flatten(fconf_lower)
    df['conf_upper'] = torch.flatten(fconf_upper)
    df['conflik_lower'] = torch.flatten(predictions_conf_lower)
    df['conflik_upper'] = torch.flatten(predictions_conf_upper)
    df['type'] = "fitted"
    # mean_const = np.repeat(torch.cat([list(m.parameters())[1].detach() for m in model.models]), len(test_x))
    # df['mean_const'] = mean_const
    print(df)
    
    # Make data frame of training input
    x_train = [(torch.stack(submodel.train_inputs).flatten()) for submodel in model.models]
    n_xtrain_cond = [len(x) for x in x_train]
    y_train = [submodel.train_targets for submodel in model.models]
    [len(y) for y in y_train]
    
    
    inputs_df =  pd.DataFrame({'peptide' : np.repeat(np.repeat(peptides2test, n_cond), n_xtrain_cond),
                            'condition' : np.repeat(np.tile(conds, n_pep), n_xtrain_cond),
                            't' : torch.cat(x_train),
                            'y': torch.cat(y_train)})
    inputs_df[inputs_df.condition != 'full model'] # remove repetitions
    inputs_df['type'] = "measured"
    df_joint = pd.concat([inputs_df, df])
    
    df_joint.to_csv(path_output + "solution_" + sm + "_" + trypticity + ".csv", sep=',', index=False, header=True)
