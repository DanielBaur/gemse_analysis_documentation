
# This Python3 library contains code focussing on the analysis of GeMSE stuff.



###############################################################
### Imports
###############################################################

import subprocess
import numpy as np
import os
import matplotlib.pyplot as plt
# include the following lines for local analysis
import getpass





###############################################################
### Generic Definitions
###############################################################


# local GeMSE analysis infrastructure installation
abspath_gemse_analysis_infrastructure = "/home/daniel/Desktop/arbeitsstuff/20200624__gemse/studies/20210105__local_gemse_analysis_framework/"
abspath_root = abspath_gemse_analysis_infrastructure +"root/root_v6.22.06.Linux-ubuntu18-x86_64-gcc7.5/root/"
pathstring_thisroot = abspath_root +"bin/thisroot.sh"
abspath_bat = abspath_gemse_analysis_infrastructure + "bat/bat/"
abspath_gemse_root_scripts = abspath_gemse_analysis_infrastructure + "gemse_root_scripts/GeMSE_ROOT_scripts/"
abspath_gemse_analysis = abspath_gemse_analysis_infrastructure + "gemse_analysis/GeMSE_analysis/"


# color
gemse_mint = "#5a8fa3" # mint green color of the GeMSE logo







###############################################################
### Helper Functions
###############################################################


# This function is used to save the 'input_data_raw_cut' structured array as a (pseudo) list file, which can further on be processed by the GeMSE analysis infrastructure
def gen_pseudo_list_file(
    pathstring_output, # pathstring according to which the pseudo list file is supposed to be saved
    header_list, # list containing the integer values of the raw data list file header
    input_data_raw_cut): # raw data (in the form of a numpy structured array) retrieved from the raw data list file via 'monxeana.get_raw_data_from_list_file()'

    # writing the file
    with open(pathstring_output, "w") as list_file:

        # writing the header of the file
        list_file.write("HEADER0:" +str(header_list[0]) +"\n")
        list_file.write("HEADER1:" +str(header_list[1]) +"\n")
        list_file.write("HEADER2:" +str(header_list[2]) +"\n")
        list_file.write("HEADER3:" +str(header_list[3]) +"\n")
        list_file.write("HEADER4:" +str(header_list[4]))

        # writing the modified, i.e. cut, data
        for i in range(len(input_data_raw_cut)):
            list_file.write("\n" +str(int(+input_data_raw_cut[i]["timestamp_ps"]/10000)) +" " +str(int(input_data_raw_cut[i]["pulse_height_adc"])) +" " +str(input_data_raw_cut[i]["flag_mca"]) +" ")

    print(f"\ngen_pseudo_list_file(): wrote '{pathstring_output}'")
    return


# This functions is used to check whether two files match line by line.
def compare_files_line_by_line(pathstring_file_a, pathstring_file_b):

    # opening both lines
    with open(pathstring_file_a, 'r') as file_a, open(pathstring_file_b, 'r') as file_b:
        ctr_lines = 0
        ctr_different_lines = 0
        # looping over both lines simultaneously utilizing 'zip'
        for line_a, line_b in zip(file_a, file_b):
            if line_a == line_b:
                continue
            else:
                print(f"lines {ctr_lines} are NOT identical: ")
                print(f"file_a: {line_a}")
                print(f"file_b: {line_b}")
                ctr_different_lines += 1
            ctr_lines += 1
    print(f"compare_files_line_by_line: {ctr_different_lines} different lines")
    return


# This function is used to extract both the base and the exponent of a float represented as a string in scientific notation (e.g., 4.77e-7).
def det_base_and_exponent_of_scientific_notation(input_float, input_precision=5):

    # formatting the input float
    if type(input_float) == float:
        input_float_string = format(input_float, "." +str(input_precision) +"e")
    elif type(input_float) == str:
        input_float_string = input_float

    # extracting base and exponent
    float_list = list(input_float_string.split("e"))
    base = float(float_list[0])
    expo = int(float_list[1])

    return base, expo


# This function is used to match the bases and exponents of a list of floats to the precision of the scientific representation of an input mean.
# The output of this function is a list of two tuples (decimal base and exponent) of all input numbers (with the mean base and exponent als the 0th element).
def match_exponents_and_precision_to_mean(
    input_mean, # input mean, either a float or a string representing a float in scientific notation
    mean_precision, # output decimal precision of the input mean (int)
    adapt_float_list): # list of numbers that has to be adapted to the input mean precision

    # formatting the input mean
    if type(input_mean) == float:
        m = format(input_mean, "." +str(mean_precision) +"e")
    elif type(input_mean) == str:
        m = format(float(input_mean), "." +str(mean_precision) +"e")

    # extracting both base and exponent for all numbers
    m_base, m_expo = det_base_and_exponent_of_scientific_notation(m)
    conv_list = [det_base_and_exponent_of_scientific_notation(adapt_float) for adapt_float in adapt_float_list]
    
    # adapting both base and exponent for all numbers
    output_list = [[format(m_base,"." +str(mean_precision) +"f"), m_expo]]
    for i in range(len(conv_list)):
        diff_expos = m_expo -conv_list[i][1]
        if diff_expos > 0:
            new_base = conv_list[i][0]/(10**(abs(diff_expos)))
        elif diff_expos <= 0:
            new_base = conv_list[i][0]*(10**(abs(diff_expos)))
        else:
            new_base = conv_list[i][0]
        #output_list.append([float(format(new_base,"." +str(mean_precision) +"f")), m_expo])
        output_list.append([format(new_base,"." +str(mean_precision) +"f"), m_expo])

    # returning a list containing all adjusted representations (i.e., a list of two tuples of the following form [<base_in_specified_precision>, <exponent_adapted_to_mean_exponent>])
    return output_list


# This function is used to convert a python isotope string (e.g., 'U238') to a latex-interpretable string (excluding the enclosing '$' characters) without requiring the rex module.
def conv_isotope_string_to_latex_syntax(isotope_string):
    amu = "".join([s for s in isotope_string if s.isdigit()])
    iso = "".join([s for s in isotope_string if s.isalpha()]).capitalize()
    latex_string = r"^{" +f"{str(int(amu))}" +r"}\mathrm{" +f"{iso}" +r"}"
    return latex_string


# This functions is used to convert a float encoded in scientific notation into a latex-interpretable string (excluding the enclosing '$' characters) without requiring the rex module
def conv_scifloat_string_to_latex_syntax(scifloat_string):
    base = float(list(scifloat_string.split("e"))[0])
    exp = int(list(scifloat_string.split("e"))[1])
    latex_string = f"{base:.2f}" +r"\cdot 10^{" +f"{exp}" +r"}"
    return latex_string





###############################################################
### GeMSE analysis infrastructure: wrapper functions
###############################################################


# This function is a wrapper for Moritz' C++ executable 'make_rootfile_list'.
def make_rootfile_list(
    input_pathstring_mca_list_file,
    input_pathstring_calibration_function,
    input_abspath_gemse_root_scripts = abspath_gemse_root_scripts):

    # executing the 'make_rootfile_list' executable
    execstring = input_abspath_gemse_root_scripts +"make_rootfile_list" +" " +input_pathstring_mca_list_file +" " +input_pathstring_calibration_function
    print(f"make_rootfile_list(): executing\n{execstring}\n")
    subprocess.call(execstring, shell=True)

    return


# This function is a wrapper for Moritz' C++ executable 'make_spectrum_list'.
def make_spectrum_list(
    input_pathstring_root_file = "XXX",
    input_abspath_gemse_root_scripts = abspath_gemse_root_scripts,
    input_time_window = [0,0]):

    # executing the 'make_spectrum_list' executable
    if input_time_window == [0,0]:
        execstring = input_abspath_gemse_root_scripts +"make_spectrum_list" +" --file " +input_pathstring_root_file +" --energy"
    else:
        execstring = input_abspath_gemse_root_scripts +"make_spectrum_list" +" --file " +input_pathstring_root_file +" --energy --t_min " +str(int(input_time_window[0])) +" --t_max " +str(int(input_time_window[1]))
    print(f"make_spectrum_list(): executing\n{execstring}\n")
    subprocess.call(execstring, shell=True)

    return


# This function is a wrapper for Moritz' C++ executable 'plot_rate'.
def plot_rate(
    input_pathstring_root_file = "XXX", # pathstring referring to the utilized root file
    input_abspath_gemse_root_scripts = abspath_gemse_root_scripts,
    input_energy_calibration = False, # flag indicating whether the energy values are supposed to be used
    input_binwidth = 1, # binwidth in seconds
    input_pulse_height_range = [0,16383]): # pulse height

    # executing the 'plot_rate' executable
    execstring = input_abspath_gemse_root_scripts +"plot_rate" +" --file " +input_pathstring_root_file
    if input_energy_calibration == True:
        execstring = execstring + " --energy"
    execstring = execstring +f" --range_min {input_pulse_height_range[0]} --range_max {input_pulse_height_range[1]}"
    execstring = execstring +f" --binwidth {input_binwidth}"
    print(f"plot_rate(): executing\n{execstring}\n")
    subprocess.call(execstring, shell=True)

    return


# This function is a wrapper for Moritz' C++ executable 'add_spectra'.
def add_spectra(
    input_pathstrings_cut_spectra,
    input_pathstring_output_spectrum,
    input_abspath_gemse_root_scripts = abspath_gemse_root_scripts):

    # executing the 'add_spectra' executable
    execstring = input_abspath_gemse_root_scripts +"add_spectra"
    for i in range(len(input_pathstrings_cut_spectra)):
        execstring = execstring +" " +input_pathstrings_cut_spectra[i]
    execstring = execstring +" " +input_pathstring_output_spectrum
    print(f"add_spectra(): executing\n{execstring}\n")
    subprocess.call(execstring, shell=True)

    return


# This function is a wrapper for Moritz' C++ executable 'GeMSE_analysis'.
def gemse_analysis(
    input_pathstring_gemse_analysis_configuration_file = "XXX",
    input_abspath_gemse_analysis = abspath_gemse_analysis):

    # executing the 'make_rootfile_list' executable
    execstring = input_abspath_gemse_analysis +"GeMSE_analysis" +" " +input_pathstring_gemse_analysis_configuration_file
    print(f"gemse_analysis(): executing\n{execstring}\n")
    subprocess.call(execstring, shell=True)

    return





###############################################################
### GeMSE analysis infrastructure: automatization
###############################################################


# This function is used to generate a `analysis_configuration_file'.
def gen_analysis_configuration_file(
    pathstring_output,
    sample_name,
    abspath_isotope_parameters_folder,
    abspath_sample_spectrum_root_file,
    abspath_background_spectrum_root_file,
    abspath_efficiency_root_file,
    abspath_resolution_root_file,
    abspath_results_folder,
    accuracy_of_mcmc = ["low", "medium", "high"][0],
    threshold_on_bayes_factor = 0.33,
    cl_for_activity_limit = 0.95,
    fract_uncert_efficency = 0.02,
    list_isotopes_to_analyze = ["U238", "Ra226", "Th228", "Ra228", "Co60", "K40", "Cs137", "Mn54", "Ti44", "Na22", "Al26"]):

    # opening and writing the output
    with open(pathstring_output, 'w+') as f:

        f.write("# sample name\n")
        f.write(f"{sample_name}\n")

        f.write("# accuracy of MCMC (low/medium/high)\n")
        f.write(f"{accuracy_of_mcmc}\n")

        f.write("# threshold on Bayes Factor for signal detection\n")
        f.write(f"{threshold_on_bayes_factor}\n")

        f.write("# CL for activity limit\n")
        f.write(f"{cl_for_activity_limit}\n")

        f.write("# fractional uncertainty of efficiencies\n")
        f.write(f"{fract_uncert_efficency}\n")

        f.write("# isotope parameters folder\n")
        f.write(f"{abspath_isotope_parameters_folder}\n")

        f.write("# sample spectrum file name\n")
        f.write(f"{abspath_sample_spectrum_root_file}\n")

        f.write("# background spectrum file name\n")
        f.write(f"{abspath_background_spectrum_root_file}\n")

        f.write("# efficiency file name\n")
        f.write(f"{abspath_efficiency_root_file}\n")

        f.write("# resolution file name\n")
        f.write(f"{abspath_resolution_root_file}\n")

        f.write("# results folder\n")
        f.write(f"{abspath_results_folder}\n")

        f.write("# isotopes to analyze\n")
        for i in range(len(list_isotopes_to_analyze)):
            f.write(f"{list_isotopes_to_analyze[i]}\n")

    return pathstring_output


# This function is used to generate a 'gemse_analysis_wiki_syntax' file from a 'gemse_analysis_summary' file.
# The idea here is to simply copy the contents of the output file into the corresponding wiki note without having to adjust every single entry manually.
def gen_analysis_results_wiki_syntax_file(pathstring_gemse_analysis_summary, pathstring_gemse_analysis_summary_wiki_syntax):

    # opening both files
    with open(pathstring_gemse_analysis_summary, 'r') as input_file, open(pathstring_gemse_analysis_summary_wiki_syntax, 'w+') as output_file:

        # looping over the contents of the input file
        flag_isotope_results = False
        flag_everything_went_fine = True
        write_string_analysis_parameters = ""
        write_string_isotopes = ""
        write_string_activity = ""
        write_string_bayes_factor = ""
        for i, line in enumerate(input_file, start=0):

            line_list = list(line.split())

            # extracting: analysis parameters
            if flag_isotope_results == False and i != 0 and "#################################" not in line and line != "\n" and not ("Isotope" in line and "Activity (Bq)" in line and "Bayes Factor" in line):
                if "sample spectrum" in line or "background spectrum" in line or "simulated efficiencies" in line or "energy resolution" in line:
                    add = line_list[0] +" " +line_list[1] +" " +"''%%" +list(line_list[2].split("/"))[-1] +"%%''"
                    print(line_list)
                else:
                    add = line[:-1]
                write_string_analysis_parameters = write_string_analysis_parameters  +add +r" \\ "

            # extracting: isotope limits/activities
            if flag_isotope_results == True:
                isotope = line_list[0]
                bayes_factor = line_list[-1]
                # case 1: limit placed
                if len(line_list) == 4 and "<" in line_list:
                    limit = line_list[2]
                    print(f"isotope: {isotope};   limit: {limit} Bq;   bayes factor: {bayes_factor}")
                    write_string_isotopes = write_string_isotopes +isotope +r" \\ "
                    write_string_activity = write_string_activity +f"< {limit}" +r" \\ "
                    write_string_bayes_factor = write_string_bayes_factor +bayes_factor +r" \\ "
                # case 2: activity determined
                elif len(line_list) == 7 and "-" in line_list and "+" in line_list:
                    activity = line_list[1]
                    activity_e_lower = line_list[3]
                    activity_e_upper = line_list[5]
                    print(f"isotope: {isotope};   activity: ({activity}-{activity_e_lower}+{activity_e_upper}) Bq;   bayes factor: {bayes_factor}")
                    write_string_isotopes = write_string_isotopes +isotope +r" \\ "
                    write_string_activity = write_string_activity +f"({activity} - {activity_e_lower} + {activity_e_upper})" +r" \\ "
                    write_string_bayes_factor = write_string_bayes_factor +bayes_factor +r" \\ "
                # case 3: exception caught
                else:
                    flag_everything_went_fine = False
                    print(f"gen_analysis_results_wiki_syntax_file(): ERROR reading line '{line_list}'")
                # writing the isotope results to the output file
                

            # checking whether the the current lines contain information on the isotope analysis
            if "Isotope" in line and "Activity (Bq)" in line and "Bayes Factor" in line:
                flag_isotope_results = True

        # printing the analysis parameters to the output file
        output_file.write(f"| GeMSE analysis | {write_string_analysis_parameters[:-4]} |||\n")
        output_file.write(f"| ::: | isotope | activity limit / measured activity [Bq] | bayes factor |\n")
        output_file.write(f"| ::: | {write_string_isotopes[:-4]} | {write_string_activity[:-4]} | {write_string_bayes_factor[:-4]} |\n\n")

        # printing the result of this function to screen
        if flag_everything_went_fine == True:
            print(f"gen_analysis_results_wiki_syntax_file(): successfully saved '{pathstring_gemse_analysis_summary_wiki_syntax}'")
        else:
            print(f"gen_analysis_results_wiki_syntax_file(): saved '{pathstring_gemse_analysis_summary_wiki_syntax}' with ERROR")

    return pathstring_gemse_analysis_summary_wiki_syntax


# This function is meant to summarize the wrapper functions for Moritz' C++ scripts defined above allowing for a all in one GeMSE analysis.
def all_in_one_gemse_analysis(
    input_pathstrings_mca_list_files, # list of pathstrings referring to the raw mca list files (all within the same folder!)
    input_time_windows, # list of time windows corresponding to the entries of 'pathstrings_mca_list_files'
    input_pathstring_calibration_function, # pathstring referring to the utilized energy calibration function
    input_pathstring_gemse_analysis_configuration_file, # pathstring referring to the analysis configuration file
    input_abspath_gemse_root_scripts = abspath_gemse_root_scripts, # abspath of Moritz' 'GeMSE_ROOT_scripts' scripts
    input_abspath_gemse_analysis = abspath_gemse_analysis): # abspath of Moritz' 'GeMSE_analysis' scripts

    ### start
    abspath_measurement_folder = input_pathstrings_mca_list_files[0][:input_pathstrings_mca_list_files[0].rfind("/")+1]
    print(f"all_in_one_gemse_analysis(): abspath_measurement_folder:\n{abspath_measurement_folder}\n")
    sepstring = "#################################################################\n"

    ### mca list file(s) ---> root file(s)
    print(sepstring +f"all_in_one_gemse_analysis(): mca list file(s) ---> root file(s)\n" +sepstring)
    for i in range(len(input_pathstrings_mca_list_files)):
        make_rootfile_list(
            input_pathstring_mca_list_file = input_pathstrings_mca_list_files[i],
            input_abspath_gemse_root_scripts = input_abspath_gemse_root_scripts,
            input_pathstring_calibration_function = input_pathstring_calibration_function)

    ### root file(s) ---> (cut) energy spectrum/spectra
    print(sepstring +f"all_in_one_gemse_analysis(): root file(s) ---> (cut) energy spectrum/spectra\n" +sepstring)
    print(input_pathstrings_mca_list_files)
    for i in range(len(input_pathstrings_mca_list_files)):
        print(input_time_windows[i])
        for j in range(len(input_time_windows[i])):
            make_spectrum_list(
                input_pathstring_root_file = input_pathstrings_mca_list_files[i] +".root",
                input_abspath_gemse_root_scripts = input_abspath_gemse_root_scripts,
                input_time_window = input_time_windows[i][j])

    ### (cut) energy spectra ---> added energy spectra
    print(sepstring +f"all_in_one_gemse_analysis(): (cut) energy spectra ---> added energy spectra\n" +sepstring)
    for i in range(len(input_pathstrings_mca_list_files)):
        # case 1: no 'time_window' was specified (i.e. '[0,0]')
        if input_time_windows[i] == [[0,0]]:
            uncut_spectrum_candidates = [abspath_measurement_folder +filename for filename in os.listdir(abspath_measurement_folder) if input_pathstrings_mca_list_files[i] +".root_spectrum_calibrated_0-" in abspath_measurement_folder +filename and filename.endswith(".root")] # the string appended when 'make_spectrum_list' is run that is indicating the live time of the spectrum is not known a priori --> accordingly one needs to search for the correct file
            if len(uncut_spectrum_candidates) == 1:
                execstring = "cp " +uncut_spectrum_candidates[0] +" " +input_pathstrings_mca_list_files[i] +".root_spectrum_calibrated_added_spectrum.root"
                print("\n", execstring, "\n")
                subprocess.call(execstring, shell=True)
            else:
                exception_string = f"all_in_one_gemse_analysis(): you specified no time cuts for file {input_pathstrings_mca_list_files[i]} and yet the following candidate files have been found:\n"
                exception_string = exception_string +''.join(["\t-->" +entry +"\n" for entry in uncut_spectrum_candidates])
                raise Exception(exception_string)
        # case 2: exactly one (non default) 'time_window' was specified
        elif len(input_time_windows[i]) == 1:
            execstring = "cp " +input_pathstrings_mca_list_files[i] +f".root_spectrum_calibrated_{str(int(input_time_windows[i][0][0]))}-{str(int(input_time_windows[i][0][1]))}s.root" +" " +input_pathstrings_mca_list_files[i] +".root_spectrum_calibrated_added_spectrum.root"
            print("\n", execstring, "\n")
            subprocess.call(execstring, shell=True)
        # case 3: multiple 'time_windows' are specified
        else:
            add_spectra(
                input_pathstrings_cut_spectra = [input_pathstrings_mca_list_files[i] +f".root_spectrum_calibrated_{str(int(itw[0]))}-{str(int(itw[1]))}s.root" for itw in input_time_windows[i]],
                input_pathstring_output_spectrum = input_pathstrings_mca_list_files[i] +".root_spectrum_calibrated_added_spectrum",
                input_abspath_gemse_root_scripts = input_abspath_gemse_root_scripts)
    # adding all added spectra together to one final spectrum that will be used for the analysis
    if len(input_pathstrings_mca_list_files) == 1:
        execstring = "cp " +input_pathstrings_mca_list_files[0] +".root_spectrum_calibrated_added_spectrum.root" +" " +abspath_measurement_folder +"final_calibrated_added_spectrum.root"
        print("\n", execstring, "\n")
        subprocess.call(execstring, shell=True)
    else:
        add_spectra(
            input_pathstrings_cut_spectra = [pathstring +".root_spectrum_calibrated_added_spectrum.root" for pathstring in input_pathstrings_mca_list_files],
            input_pathstring_output_spectrum = abspath_measurement_folder +"final_calibrated_added_spectrum",
            input_abspath_gemse_root_scripts = input_abspath_gemse_root_scripts)

    ### bayesian analysis
    print(sepstring +f"all_in_one_gemse_analysis(): bayesian analysis\n" +sepstring)
    # printing the analysis settings
    print(f"all_in_one_gemse_analysis(): analysis settings")
    sample_name = ""
    results_folder = ""
    with open(input_pathstring_gemse_analysis_configuration_file) as analysis_settings_file:
        for line in analysis_settings_file:
            print("\t", line[:-1])
            if sample_name == "active":
                sample_name = line[:-1]
            if "# sample name" in line:
                sample_name = "active"
            if results_folder == "active":
                results_folder = line[:-1]
            if "# results folder" in line:
                results_folder = "active"
    print("")
    print(f"all_in_one_gemse_analysis(): results_folder='{results_folder}'\n")
    print(f"all_in_one_gemse_analysis(): sample_name='{sample_name}'\n")
    print("")
    # running the analysis
    gemse_analysis(
        input_pathstring_gemse_analysis_configuration_file = input_pathstring_gemse_analysis_configuration_file,
        input_abspath_gemse_analysis = input_abspath_gemse_analysis)

    ### aftermath
    # printing the analysis results
    print(f"all_in_one_gemse_analysis(): analysis results")
    with open(results_folder +sample_name +"_activities_summary.txt") as analysis_results_file:
        for line in analysis_results_file:
            print("\t", line[:-1])
    print("")

    ### end
    return



    

###############################################################
### software-based veto investigation
###############################################################


# This is the dtype used for raw data extracted from MC2Analyzer.
timestamp_data_mc2_dtype = np.dtype([
    ("timestamp_10ns", np.uint64), # timestamp in 10ns
    ("pulse_height_adc", np.int64), # max adc channel is ~16000, np.int16 ranges from -32768 to 32767
    ("extra", np.int32), # 
    ("validity", np.unicode_, 16), # 
])


# This function is used to load the list file generated by the MCA.
def get_timestamp_data_as_ndarray(pathstring_data = ""):
    timestamp_data_tuplelist = []
    fname = "get_timestamp_data_as_ndarray"
    with open(pathstring_data) as input_file:
        for line in input_file:
            if not line.startswith("HEADER"):
                line_list = list(line.split())
                try:
                    timestamp_10ns = np.uint64(line_list[0])
                    pulse_height_adc = np.int64(line_list[1])
                    extra = np.int32(line_list[2])
                    timestamp_data_tuplelist.append((
                        timestamp_10ns,
                        pulse_height_adc,
                        extra,
                        "valid"))
                except:
                    print(fname, ": exception: ", line_list)
    return np.array(timestamp_data_tuplelist, timestamp_data_mc2_dtype)


# This function is used to compare both the signal and veto file generated by the MCA.
# According to the veto file the entries in the signal file are then either marked as "vetoed" or "not_vetoed".
def get_veto_information(
    input_signal_file,
    pathstring_vetodata,
    timingoffset = 10, # in us
    vetowindow = 10):# in us

    # initial definitions
    signal_file = input_signal_file.copy()
    exception_list = []
    fname = "get_veto_information"
    t_i = datetime.datetime.now()
    print(fname, " : started: ", t_i)
    print(fname, " : input veto file: ", pathstring_vetodata)
    l = len(signal_file)
    o = timingoffset*100 # the timestamp recorded by the MCA corresponds to clock cycles, i.e. 10ns
    v = vetowindow*100 # accordingly one must convert us to 10ns
    j = 0 # index of the current signal file entry to be checked

    # accessing and looping over the veto file line by line (note that therefore the file does not have to be loaded into the RAM in its entirety)
    with open(pathstring_vetodata) as input_file:
        for line in input_file:
            if not line.startswith("HEADER"):
                line_list = list(line.split())
                try:
                    timestamp_10ns = np.uint64(line_list[0])
                    pulse_height_adc = np.int64(line_list[1])
                    # go to the next veto entry if the current signal entry timestamp is larger than the current veto entry timestamp
                    if signal_file[j]["timestamp_10ns"] > timestamp_10ns +o +v:
                        continue
                    # if the veto entry timestamp is larger than the current signal entry, check whether the signal entry timestamp lies within the interval [timestamp_10ns +o, timestamp_10ns +o +w] and therefore needs to be vetoed, otherwiese bring up the next signal entry until their timestamp is greater than timestamp_10ns +o +w (and one would again have to skip lines until a smaller signal entry timestamp is once again found)
                    else:
                        while signal_file[j]["timestamp_10ns"] <= timestamp_10ns +o +v and j<l:
                            if signal_file[j]["timestamp_10ns"] > timestamp_10ns +o:
                                if signal_file[j]["validity"] == "cut":
                                    signal_file[j]["validity"] = "cut_and_vetoed"
                                    print(fname, " : vetoed already cut entry: ", signal_file[j])
                                else:
                                    signal_file[j]["validity"] = "vetoed"
                                    print(fname, " : vetoed: ", signal_file[j])
                                j +=1
                            else:
                                j +=1
                except:
                    exception_list.append(line_list)
            if j % 10000 == 0:
                print(fname, ": ", f"checked {j} of {l} entries, already running for {datetime.datetime.now()-t_i}" )
    t_f = datetime.datetime.now()
    print(fname, ": ", "encountered exceptions:")
    for i in range(len(exception_list)):
        print(exception_list[i])
    print(fname, ": ", "finished: ", t_f)
    print(fname, ": ", "duration: ", t_f-t_i)

    return signal_file


# This function is used to cut data from the input signal file.
def get_cut_information(input_signal_file):

    signal_file = input_signal_file.copy()
    for i in range(len(signal_file)):
        if (signal_file[i]["extra"]>8 and signal_file[i]["extra"]<16) or signal_file[i]["pulse_height_adc"]>17000 or signal_file[i]["pulse_height_adc"]<0:
            if signal_file[i]["validity"] == "vetoed":
                signal_file[i]["validity"] = "cut_and_vetoed"
            else:
                signal_file[i]["validity"] = "cut"
    return signal_file


# This function is used to display an overview of the signal file
def display_signal_file_properties(signal_file):
    l = len(signal_file)
    l_valid = len(signal_file[(signal_file["validity"]=="valid")])
    l_veto = len(signal_file[(signal_file["validity"]=="vetoed")|(signal_file["validity"]=="cut_and_vetoed")])
    l_cut = len(signal_file[(signal_file["validity"]=="cut")|(signal_file["validity"]=="cut_and_vetoed")])
    l_cutveto = len(signal_file[(signal_file["validity"]=="cut_and_vetoed")])
    print(f"signal file properties:")
    print(f"    number of events: {l}")
    print(f"    valid: {l_valid} ({l_valid/l*100}%)")
    print(f"    vetoed: {l_veto} ({l_veto/l*100}%)")
    print(f"    cut: {l_cut} ({l_cut/l*100}%)")
    print(f"    cutveto: {l_cutveto} ({l_cutveto/l*100}%)")
    return


#if [True, False][1]:
#
#
#    ### definitions
#    # small files
#    filestring_raw0 = "/media/daniel/intenso_256gb_stick/2020-07-07__ptfe_02__ek_clean/2020-08-04__ptfe_02_ek_clean.txt_ch000-20200806_1423_copy.txt"
#    filestring_raw1 = "/media/daniel/intenso_256gb_stick/2020-07-07__ptfe_02__ek_clean/2020-08-04__ptfe_02_ek_clean.txt_ch001-20200806_1423_copy.txt"
#    # big files
#    filestring_raw0 = "/media/daniel/intenso_256gb_stick/2020-07-07__ptfe_02__ek_clean/2020-08-04__ptfe_02_ek_clean.txt_ch000-20200815_1649_copy.txt"
#    filestring_raw1 = "/media/daniel/intenso_256gb_stick/2020-07-07__ptfe_02__ek_clean/2020-08-04__ptfe_02_ek_clean.txt_ch001-20200815_1649_copy.txt"
#
#
#    # loading raw data
#    ch0_raw = get_timestamp_data_as_ndarray(filestring_raw0) # len(ch0_raw) = 77455
#
#
#    # veto
#    ch0_veto = get_veto_information(# len(ch1_raw) = 121307796
#        input_signal_file = ch0_raw,
#        pathstring_vetodata = filestring_raw1,
#        timingoffset = -10, # in us
#        vetowindow = 10)# in us
#
#
#    # cuts (depending on "extra")
#    ch0_veto_cut = get_cut_information(input_signal_file=ch0_veto)
#
#
#    # displaying cut and veto results
#    display_signal_file_properties(signal_file=ch0_veto_cut)





###############################################################
### PTFEsc-specific analysis stuff
###############################################################


if getpass.getuser() == "daniel":

    import uproot
    import sys
    sys.path.append("/home/daniel/Desktop/arbeitsstuff/20180705__monxe/monxe_software/miscfig/")
    sys.path.append("/home/daniel/Desktop/arbeitsstuff/20180705__monxe/monxe_software/monxeana/") 
    import monxeana
    import Miscellaneous_Figures as miscfig
    import json


    # This function is used to .
    def gemse_analysis_aftermath(
        input_filenames, # list containing the utilized mca list files
        input_time_windows, # list containing the utilized time windows for the respective mca list files
        input_sample_mass, # mass of the examined sample in kg
        input_pathstring_calibration_function, # pathstring referring to the utilized calibration function
        input_pathstring_gemse_analysis_summary, # pathstring referring to the analysis output summary file
        input_pathstring_added_root_spectrum, # added spectrum root file
        input_pathstring_wiki_syntax_output, # pathstring referring to the ouput wiki syntax file (which is supposed to simply be copied into the PTFEsc wiki note)
        input_pathstring_json_output, # pathstring referring to the output .json file
        input_pathstrings_spectrum_plot, # pathstrings referring to the output commented spectrum plot
        flag_config = ["default"][0],
        input_ylim = ""): # keywords passed on to the 'ax1.set_ylim()' function call

        """
        This function is used to provide a all-in-one function call automatically generating an elaborate summarizing output of the GeMSE analysis of a specific sample.
        I.e., wiki syntax output, commented spectrum plot.
        """
        
        ### storing the information from the analysis summary file in a dictionary
        analysis_dictionary = {
            "mca_list_files" : input_filenames,
            "mca_list_files_time_windows" : input_time_windows,
            "sample_mass_kg" : input_sample_mass,
            "calibration_function" : list(input_pathstring_calibration_function.split("/"))[-1],
            "isotope_data" : {}}
        with open(input_pathstring_gemse_analysis_summary, 'r') as gemse_analysis_summary_file:
            flag_isotope_results = False
            for i, line in enumerate(gemse_analysis_summary_file, start=1):
                line_list = list(line.split())
                # extracting parameters
                if i==2:
                    analysis_dictionary.update({"datetimestamp" : line_list[0] +" " +line_list[1]})
                elif "sample name" in line:
                    analysis_dictionary.update({"sample_id" : list(line_list[2].split("/"))[-1]})
                elif "sample spectrum" in line:
                    analysis_dictionary.update({"sample_spectrum" : list(line_list[2].split("/"))[-1]})
                elif "background spectrum" in line:
                    analysis_dictionary.update({"background_spectrum" : list(line_list[2].split("/"))[-1]})
                elif "simulated efficiencies" in line:
                    analysis_dictionary.update({"simulated_efficiencies" : list(line_list[2].split("/"))[-1]})
                elif "fractional uncertainty efficiencies" in line:
                    analysis_dictionary.update({"fractional_uncertainty_efficiencies" : line_list[3]})
                elif "energy resolution" in line:
                    analysis_dictionary.update({"energy_resolution" : list(line_list[2].split("/"))[-1]})
                elif "measurement time sample" in line:
                    analysis_dictionary.update({"measurement_time_s" : line_list[3]})
                    analysis_dictionary.update({"measurement_time_d" : float(line_list[3]) / (60*60*24)})
                elif "measurement time background" in line:
                    analysis_dictionary.update({"measurement_time_background_sec" : line_list[3]})
                elif "BF threshold for signal" in line:
                    analysis_dictionary.update({"bf_threshold_for_signal" : line_list[4]})
                elif "CL for activity limit" in line:
                    analysis_dictionary.update({"cl_for_activity_limit" : line_list[4]})
                # extracting isotope limits/activities
                elif flag_isotope_results == True:
                    isotope = line_list[0]
                    analysis_dictionary["isotope_data"].update({
                        isotope : {
                            "bayes_factor" : line_list[-1],
                            "upper_limit_bq" : "",
                            "activity_bq" : "",
                            "activity_bq_lower" : "",
                            "activity_bq_upper" : "",
                            "upper_limit_bq_per_kg" : "",
                            "activity_bq_per_kg" : "",
                            "activity_bq_lower_per_kg" : "",
                            "activity_bq_upper_per_kg" : ""}})
                    if len(line_list) == 4 and "<" in line_list:
                        analysis_dictionary["isotope_data"][isotope]["upper_limit_bq"] = float(line_list[2])
                        analysis_dictionary["isotope_data"][isotope]["upper_limit_bq_per_kg"] = float(line_list[2]) / analysis_dictionary["sample_mass_kg"]
                    elif len(line_list) == 7 and "-" in line_list and "+" in line_list:
                        analysis_dictionary["isotope_data"][isotope]["activity_bq"] = line_list[1]
                        analysis_dictionary["isotope_data"][isotope]["activity_bq_lower"] = line_list[3]
                        analysis_dictionary["isotope_data"][isotope]["activity_bq_upper"] = line_list[5]
                        analysis_dictionary["isotope_data"][isotope]["activity_bq_per_kg"] = float(line_list[1]) / analysis_dictionary["sample_mass_kg"]
                        analysis_dictionary["isotope_data"][isotope]["activity_bq_lower_per_kg"] = float(line_list[3]) / analysis_dictionary["sample_mass_kg"]
                        analysis_dictionary["isotope_data"][isotope]["activity_bq_upper_per_kg"] = float(line_list[5]) / analysis_dictionary["sample_mass_kg"]
                    else:
                        raise Exception(f"something went wrong: {line_list}")
                elif "Isotope" in line and "Activity (Bq)" in line and "Bayes Factor" in line:
                    flag_isotope_results = True
                else:
                    continue

        ### saving the analysis results dictionary as a .json file
        with open(input_pathstring_json_output, "w") as json_output_file:
            json.dump(analysis_dictionary, json_output_file, indent=4)

        ### generating the wiki syntax output file
        with open(input_pathstring_wiki_syntax_output, 'w+') as output_file:
            # analysis parameters
            analysis_parameters_list = [
                "sample ID:   " +"''" +analysis_dictionary['sample_id'] +"''",
                "measurement files:   " +"''" +r"'', ''".join(input_filenames) +"''",
                "analysis date:   " +analysis_dictionary['datetimestamp'],
                "energy resolution:   " +"''" +analysis_dictionary['energy_resolution'] +"''",
                "background:   " +"''" +analysis_dictionary['background_spectrum'] +"''",
                "efficiencies:   " +"''" +analysis_dictionary['simulated_efficiencies'] +"''",
                "calibration:   " +"''" +list(input_pathstring_calibration_function.split("/"))[-1] +"''",
                "fractional uncertainty efficiencies:   " +analysis_dictionary['fractional_uncertainty_efficiencies'],
                "BF threshold for signal:   " +analysis_dictionary['bf_threshold_for_signal'],
                "CL for activity limit:   " +analysis_dictionary['cl_for_activity_limit'],
            ]
            write_string_analysis_parameters = r" \\ ".join(analysis_parameters_list)
            # isotope parameters
            write_string_isotopes = ""
            write_string_activity = ""
            write_string_bayes_factor = ""
            for key in analysis_dictionary["isotope_data"].keys():
                write_string_isotopes += key +r" \\ "
                if analysis_dictionary['isotope_data'][key]['upper_limit_bq'] != "":
                    write_string_activity += f"< {analysis_dictionary['isotope_data'][key]['upper_limit_bq']}" +r" \\ "
                else:
                    write_string_activity += f"({analysis_dictionary['isotope_data'][key]['activity_bq']} - {analysis_dictionary['isotope_data'][key]['activity_bq_lower']} + {analysis_dictionary['isotope_data'][key]['activity_bq_upper']})" +r" \\ "
                write_string_bayes_factor += f"{analysis_dictionary['isotope_data'][key]['bayes_factor']}" +r" \\ "
            # printing to the output file
            output_file.write(f"| GeMSE analysis | {write_string_analysis_parameters} |||\n")
            output_file.write(f"| ::: | isotope | activity limit / measured activity [Bq] | bayes factor |\n")
            output_file.write(f"| ::: | {write_string_isotopes[:-4]} | {write_string_activity[:-4]} | {write_string_bayes_factor[:-4]} |\n\n")
        print(f"gemse_analysis_aftermath(): saved {input_pathstring_wiki_syntax_output}")

        ### generating the commented spectrum output plot
        for flag_plot in ["plain","commented_summary"]:#, "commented_internal"]:
            # extracting the data from the added spectrum root file
            added_root_spectrum = uproot.open(input_pathstring_added_root_spectrum)
            hist = added_root_spectrum["hist"]
            bin_edges = hist.axis().edges() # aequidistant engergy bin edges
            bin_centers = [bin_edges[i] +0.5*(bin_edges[i+1]-bin_edges[i]) for i in range(len(bin_edges)-1)]
            counts = list(hist.values()) # number of counts per energy bin
            counts_errors = list(hist.errors()) # 
            # figure formatting
            fig, ax1 = plt.subplots(figsize=miscfig.image_format_dict["talk"]["figsize"], dpi=150)
            #y_lim = [0, 1.1*(max(counts) +max(counts_errors))]
            x_lim = [bin_edges[0], bin_edges[-1]]
            ax1.set_xlim(x_lim)
            ax1.set_yscale('log')
            if input_ylim != "":
                ax1.set_ylim(input_ylim)
            ax1.yaxis.set_ticklabels([], minor=True)
            ax1.set_xlabel("energy deposition / $\mathrm{keV}$")
            binwidth = float(bin_centers[2]-bin_centers[1])
            ax1.set_ylabel("entries per " +f"${binwidth:.1f}" +r"\,\mathrm{keV}$")
            # plotting the stepized histogram
            bin_centers, counts, counts_errors_lower, counts_errors_upper, bin_centers_mod, counts_mod = monxeana.stepize_histogram_data(
                bincenters = bin_centers,
                counts = counts,
                counts_errors_lower = counts_errors,
                counts_errors_upper = counts_errors,
                flag_addfirstandlaststep = True)
            plt.plot(
                bin_centers_mod,
                counts_mod,
                linewidth = 0.2,
                color = "black",
                linestyle='-',
                zorder=1,
                label="jfk")
            plt.fill_between(
                bin_centers,
                counts-counts_errors_lower,
                counts+counts_errors_upper,
                color = gemse_mint,
                alpha = 1,
                zorder = 0,
                interpolate = True)
            # annotations
            if flag_plot == "commented_summary":
                comment_list_sample = [r"\texttt{" +analysis_dictionary['sample_id'].replace("_","\_") +r"} ($" +f"{analysis_dictionary['sample_mass_kg']:.1f}" +r"\,\mathrm{kg},\," +f"{analysis_dictionary['measurement_time_d']:.1f}" +r"\,\mathrm{d}" +"$)"]
                comment_list_files = [r"   \texttt{" +f.replace("_","\_") +r"}" for i,f in enumerate(input_filenames)]
                comment_list_results = [
                    r"$" +conv_isotope_string_to_latex_syntax(key) +r"$: $<" +conv_scifloat_string_to_latex_syntax(format((float(analysis_dictionary['isotope_data'][key]["upper_limit_bq"])/input_sample_mass), ".2e")) +r"\,\mathrm{Bq/kg}$" 
                    if analysis_dictionary['isotope_data'][key]["upper_limit_bq"] != "" 
#                    else "" 
                    else r"$" +conv_isotope_string_to_latex_syntax(key) +r"$: $(" +f"{match_exponents_and_precision_to_mean(analysis_dictionary['isotope_data'][key]['activity_bq_per_kg'], 2, [analysis_dictionary['isotope_data'][key]['activity_bq_upper_per_kg'], analysis_dictionary['isotope_data'][key]['activity_bq_lower_per_kg']])[0][0]}" +r"^{+" +f"{match_exponents_and_precision_to_mean(analysis_dictionary['isotope_data'][key]['activity_bq_per_kg'], 2, [analysis_dictionary['isotope_data'][key]['activity_bq_upper_per_kg'], analysis_dictionary['isotope_data'][key]['activity_bq_lower_per_kg']])[1][0]}" +r"}" +r"_{-" +f"{match_exponents_and_precision_to_mean(analysis_dictionary['isotope_data'][key]['activity_bq_per_kg'], 2, [analysis_dictionary['isotope_data'][key]['activity_bq_upper_per_kg'], analysis_dictionary['isotope_data'][key]['activity_bq_lower_per_kg']])[2][0]}" +r"})\cdot 10^{" +f"{match_exponents_and_precision_to_mean(analysis_dictionary['isotope_data'][key]['activity_bq_per_kg'], 2, [analysis_dictionary['isotope_data'][key]['activity_bq_upper_per_kg'], analysis_dictionary['isotope_data'][key]['activity_bq_lower_per_kg']])[0][1]}" +"}" +r"\,\mathrm{Bq/kg}$" 
                    for i,key in enumerate(analysis_dictionary['isotope_data'].keys())]
#mean base: match_exponents_and_precision_to_mean(analysis_dictionary['isotope_data'][key]['activity_bq_per_kg'], 2, [analysis_dictionary['isotope_data'][key]['activity_bq_upper_per_kg'], analysis_dictionary['isotope_data'][key]['activity_bq_lower_per_kg']])[0][0]
#mean expo: match_exponents_and_precision_to_mean(analysis_dictionary['isotope_data'][key]['activity_bq_per_kg'], 2, [analysis_dictionary['isotope_data'][key]['activity_bq_upper_per_kg'], analysis_dictionary['isotope_data'][key]['activity_bq_lower_per_kg']])[0][1]
#upper base: match_exponents_and_precision_to_mean(analysis_dictionary['isotope_data'][key]['activity_bq_per_kg'], 2, [analysis_dictionary['isotope_data'][key]['activity_bq_upper_per_kg'], analysis_dictionary['isotope_data'][key]['activity_bq_lower_per_kg']])[1][0]
#lower base: match_exponents_and_precision_to_mean(analysis_dictionary['isotope_data'][key]['activity_bq_per_kg'], 2, [analysis_dictionary['isotope_data'][key]['activity_bq_upper_per_kg'], analysis_dictionary['isotope_data'][key]['activity_bq_lower_per_kg']])[2][0]
                monxeana.annotate_comments(
                    comment_ax = ax1,
                    comment_list = comment_list_sample,
                    comment_textpos = [0.025, 0.930],
                    comment_textcolor = "black",
                    comment_linesep = 0.1,
                    comment_fontsize = 11)
                monxeana.annotate_comments(
                    comment_ax = ax1,
                    comment_list = comment_list_results,
                    comment_textpos = [0.970, 0.930],
                    comment_textcolor = "black",
                    comment_linesep = 0.083,
                    comment_fontsize = 9)
            elif flag_plot == "commented_internal":
                comment_list_sample = [r"\texttt{" +input_sample_id.replace("_","\_") +r"}"]
                comment_list_results = ["yey"]
            # saving the output plot
            for i in input_pathstrings_spectrum_plot:
                if i != "":
                    savepathstring = i[:-4] +"__" +flag_plot +i[-4:]
                    fig.savefig(savepathstring)
                    print(f"gemse_analysis_aftermath(): saved {savepathstring}")

        return input_pathstring_json_output


    # This function is used to nicely print the analysis results stored in the .json referred to by 'pathstring_analysis_results_dictionary'
    def print_analysis_results_nicely(pathstring_analysis_results_json_pathstring):

        # loading the .json summary file
        with open(pathstring_analysis_results_json_pathstring, "r") as json_input_file:
            analysis_results_json_file = json.load(json_input_file)
        print(f"\n\n\n\n")
        print(100*"#")
        print(f"### print_analysis_results_nicely: summarized results of {analysis_results_json_file['sample_id']}")
        print(100*"#")
        print(f"\n\n")

        # printing input
        for key in analysis_results_json_file.keys():
            if key not in ["isotope_data"]:
                print(f"{key}:")
                print(f"\t{analysis_results_json_file[key]}\n")

        # printing isotope data
        print(f"isotope_data:")
        for ke in analysis_results_json_file["isotope_data"].keys():
            print(f"\t{ke}:")
            for k in analysis_results_json_file["isotope_data"][ke].keys():
                print(f"\t\t{k}: {analysis_results_json_file['isotope_data'][ke][k]}")

        return
  


