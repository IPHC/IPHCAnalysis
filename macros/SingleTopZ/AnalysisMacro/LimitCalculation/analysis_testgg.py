# for model building:
def get_model(signalname):

    # Read in and build the model automatically from the histograms in the root file. 
    # This model will contain all shape uncertainties given according to the templates
    # which also includes rate changes according to the alternate shapes.
    # For more info about this model and naming conventuion, see documentation
    # of build_model_from_rootfile.
    model = build_model_from_rootfile('NewFileToBeUsedForThetaWithAutoNamingConvention_allpoints.root')

    # If the prediction histogram is zero, but data is non-zero, teh negative log-likelihood
    # is infinity which causes problems for some methods. Therefore, we set all histogram
    # bin entries to a small, but positive value:
    model.fill_histogram_zerobins()

    # define what the signal processes are. All other processes are assumed to make up the 
    # 'background-only' model.
    model.set_signal_processes('FCNC_'+signalname+'*')

    # Add some lognormal rate uncertainties. The first parameter is the name of the
    # uncertainty (which will also be the name of the nuisance parameter), the second
    # is the 'effect' as a fraction, the third one is the process name. The fourth parameter
    # is optional and denotes the channl. The default '*' means that the uncertainty applies
    # to all channels in the same way.
    # Note that you can use the same name for a systematic here as for a shape
    # systematic. In this case, the same parameter will be used; shape and rate changes 
    # will be 100% correlated.
    model.add_lognormal_uncertainty('Zjets_rate', math.log(1.5), 'DataZjets')
    model.add_lognormal_uncertainty('WZ_rate', math.log(1.2), 'WZ')
    
    model.add_lognormal_uncertainty('TTbar_rate', math.log(1.1), 'TTbar')
    model.add_lognormal_uncertainty('ZZ_rate', math.log(1.3), 'ZZ')
    
    model.add_lognormal_uncertainty('norm1', math.log(1.3), 'FCNC_zut05')
    model.add_lognormal_uncertainty('norm2', math.log(1.3), 'FCNC_zut10')
    model.add_lognormal_uncertainty('norm3', math.log(1.3), 'FCNC_zut20')
    model.add_lognormal_uncertainty('norm4', math.log(1.3), 'FCNC_zut30')
    model.add_lognormal_uncertainty('norm5', math.log(1.3), 'FCNC_zut50')
    model.add_lognormal_uncertainty('norm5', math.log(1.3), 'FCNC_zut74')

	
	
    # the qcd model is from data, so do not apply a lumi uncertainty on that:
    for p in model.processes:
        if p == 'WZ': continue
        if p == 'DataZjets': continue
        model.add_lognormal_uncertainty('lumi',    math.log(1.025), p)
        model.add_lognormal_uncertainty('Lept',    math.log(1.01) , p)
        model.add_lognormal_uncertainty('Trig',    math.log(1.05) , p)
        model.add_lognormal_uncertainty('pu',      math.log(1.01) , p)
        model.add_lognormal_uncertainty('BTag',    math.log(1.01) , p)
        #model.add_lognormal_uncertainty('Brs',     math.log(1.05) , p)
    # Specifying all uncertainties manually can be error-prone. You can also execute
    # a automatically generated file using python's execfile here
    # which contains these statements, or read in a text file, etc. Remember: this is a
    # python script, so use this power!
    return model

# -------------- TO CHANGE BY THE USER
signalname = 'zut'
# -------------- TO CHANGE BY THE USER
model = get_model(signalname)


# first, it is a good idea to generate a summary report to make sure everything has worked
# as expected. The summary will generate quite some information which should it make easy to spot
# errors like typos in the name of uncertainties, missing shape uncertaintie, etc.
model_summary(model)




# 2.b. CLs limits
# calculate cls limit plots. The interface is very similar to bayesian_limits. However, there are a few
# more options such as the definition of the test statistic which is usually a likelihood ratio but can differ in
# which parameters are minimized and which constraints / ranges are applied during minimization.
# Here, we stay with the default which fixes beta_signal=0
# for the background only hypothesis and lets it float freely for the signal+background hypothesis.
# See cls_limits documentation for more options.

# la ligne ci dessous ne fonctionne pas=> changer directement cls_limit.py

# ------------------- TO CHANGE BY THE USER
xsections = ['05','10','20','30','50']
# ------------------- TO CHANGE BY THE USER

xsections2 = []
for item in xsections :
    tmp=[]
    tmp.append('FCNC_'+signalname+item)
    xsections2.append(tmp)
print "Signal processes : " + str(xsections2)

plot_exp, plot_obs = cls_limits(model, signal_processes = xsections2)

#plot_exp, plot_obs = cls_limits(model, ts_column = ['lr'], signal_processes = [['FCNC_zut']])
#plot_exp, plot_obs = cls_limits(model, signal_processes = [['FCNC_zut05'],['FCNC_zut10'],['FCNC_zut20'],['FCNC_zut30'],['FCNC_zut50'],['FCNC_zut74']])
#plot_exp, plot_obs = cls_limits(model, signal_processes = [['FCNC_zut05'],['FCNC_zut10'],['FCNC_zut20']])
#plot_exp, plot_obs = cls_limits(model, signal_processes = [['FCNC_zut05']])
#plot_exp, plot_obs = cls_limits(model, signal_processes = [['FCNC_zut20']])
#plot_exp, plot_obs = cls_limits(model, ts_column = ['lhclike'], signal_processes = [['FCNC_zut']])

# as for the bayesian limits: write the result to a text file
plot_exp.write_txt('cls_limits_expected.txt')
plot_obs.write_txt('cls_limits_observed.txt')

