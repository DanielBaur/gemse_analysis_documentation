
# This Python3 library contains code focussing on the analysis of GeMSE stuff.



###############################################################
### Imports
###############################################################

import subprocess
import numpy as np
import os





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


