'''This code is meant to find the equivalent resistance of a glass sample, studied using an EIS machine.
This is meant to locate the resistance value in the dip after the semicircle of a Nyquist plot. As such,
the first local minimum from the right is located. This may detect faulty data if there is no tail recorded,
if there is no circle recorded, or if there is a random fluctuation in the data.'''

'''Note that in the local minimum locator there is a minimum to how far down in frequency one must go before
counting a local minimum, to avoid issues with stray points. This may be enabled or disabled. Also, you may
adjust which value flags for violations and which ones are prohibited entirely.'''

'''ENSURE THAT you adjust the number of data sets to take in, both the upper and lower bounds.'''


import math

#__________CHANGE BEFORE USE______(Also Change File Name Format Below)__________________________________
n_data_points_per_set = 49 # creates a variable for number of data points per data set. Initially, this was just 49. not all instances of 49 may have been changed to this variable
n_data_sets_total = 312 # creates a variable for the total number of data sets we have (including the ones you want to exclude from the beginning)
n_data_sets_exclude_at_start = 1 # creates a variable for the number of data sets we want to exclude at the beginning
ramp_data_sets = [1] #INDEXED FROM ONE creates a variable for the data sets from temperature ramps, which have nyquist plots at varying temperature values, format [1, 2, 4...]
ramp_data_sets_imaginary_column = 3 #creates a variable for the column number INDEXED FROM 0 for imaginary values
ramp_data_sets_real_column = 2 #creates a variable for the column number INDEXED FROM 0 for real values
ramp_data_sets_temp_column = 1 #creates a variable for the column number INDEXED FROM 0 for temperature values
sample_thickness_mm = 2.00 #creates a variable for the thickness of the sample in mm
sample_electrode_area_mm2 = ((4.46 / 2) ** 2 ) * 3.14159 #creates a variable for the area of the electrodes in mm^2
intended_flagged_value = 10 # sets the minimum number from the left to be considered problematic, normal indexing, from 0
intended_forbidden_value_right = 30 # sets the maximum number from the left to be allowed, inclusive, indexed from zero
intended_forbidden_value_left = 5 # sets the minimum number from the left to be allowed, inclusive, indexed from zero
low_separation = 2 #creates minimum separation of low values, on first glance I recommend 3 or 4 on first glance
heat_number = "first" #which heat is this? format as a string like "first", "second", "fourth"
heat_temperature = 484 #what temperature was this heated at? put it in celcius
sample_name = ".25 Li2O .75 B2O3" #what is the name of the sample? put it in the format of .45 Na2O .55 SiO2 as a string
min_rise_dist = 1 #minimum number of points that must be above and to the left of the local minimum, when scanning from the right
which_to_nyquist_plot = list(range(1, n_data_sets_total, 20)) #list(range(5, 433, 20))#[1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90]  #INDEXED FROM ONE creates a list of data sets by file name at which to take nyquist plots. syntax of range is start, stop, step. may also manually input list
range_to_lin_fit_cond = [[4550, 5590]]#[[2820, 3480], [3600, 4800]] #creates a list of pairs of inclusive lower and upper bounds of seconds to fit resistivity and conductivity to as a line. also, will print nyquist plots of all of those points. Formatted [[min1, max1], [min2, max2], ...]
plot_ranges_overplot = []#[[10, 20, 30, 40], [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280]] #creates a list of lists of indices (index from one in x of 2501 format) to do nyquist plots on, overplotting them together, of the form [[1, 2, 3, 4, 20], [101, 102, 105, 110, 167]]
imaginary_column = 3 # column for the imaginary impedance values in the CSV, indexed from 0
real_column = 2 # column for the real impedance values in the CSV, indexed from 0
freq_column = 0 #column for the frequency values NEW 3-24-25
#____________Fill out the which heat, which temp, what sample___________________________________________
#check out the TODOs to find other things you may need to edit
#also may need tto edit the column location of data


#____________________________error determiner_______, frequency resolution on the spec sheet
# go error region by error region
def uncertainty_re_im(re, neg_im, freq): #re is the z_real value, neg_im is negative imaginary z value, f_id is the index of the frequency, beginning at 0
    total_imp = math.sqrt(re**2 + neg_im**2) #impedance in ohms
    log_imp = math.log10(total_imp) # log of impedance in ohms
    log_freq = math.log10(freq) #log of frequency in hertz
    z_unc = 1000 #uncertainty we assign for being outside all ranges in percent
    ang_unc = 180 #uncertianty in argument in degrees
    #add the points and ranges just for each uncertainty zone
    all_zone_freqs = [[0.0000025, 0.35, 100, 33000, 1000000, 20000000],
                      [0.0000025, 0.33, 100, 33000, 1000000, 10000000],
                      [0.0000025, 0.29, 100, 33000, 1000000, 4600000],
                      [0.0000025, 0.047, 0.34, 10, 24000, 33000, 2000000],
                      [0.000205, 0.047, 0.34, 10, 600000, 2000000]
                      ] # define a list of lists of endpoints of frequency ranges, left to right in Hz
    all_zone_upper_imps = [[400000000000000, 100000000000000, 2000000000000, 4000000000, 150000000, 10000],
                           [40000000000000, 10000000000000, 200000000000, 500000000, 20000000, 10000],
                           [3500000000000, 1000000000000, 17000000000, 43000000, 1400000, 7000],
                           [1000000000, 1000000000, 130000000, 4400000, 1900, 1900, 1900],
                           [200000000000, 120000000000, 100000000000, 10000000000, 107000, 10000]
                           ] # define a list of lists of upper limits for impedance at endpoints of frequency ranges, left to right in Ohms
    all_zone_lower_imps = [[0.01, 0.01, 0.01, 0.01, 0.3, 6],
                           [0.1, 0.1, 0.1, 0.1, 3.2, 30],
                           [1, 1, 1, 1, 35, 150],
                           [10, 10, 10, 10, 10, 10, 600],
                           [1000000000, 1000000000, 130000000, 4400000, 1900, 1900]
                           ] # define a list of lists of lower limits for impedance at endpoints of frequency ranges, left to right in Ohms
    z_unc_of_zone = [100, 10, 1, 0.2, 0.1] #define a list of uncertainties in impedance magnitude in percent of the uncertainty ranges
    ang_unc_of_zone = [180, 6, 0.6, 0.06, 0.06] #define a list of uncertainties in impedance argument in degrees of the uncertainty ranges
    for zone_idx in range(len(all_zone_freqs)):
        zone_freqs = all_zone_freqs[zone_idx]
        zone_upper_imps = all_zone_upper_imps[zone_idx]
        zone_lower_imps = all_zone_lower_imps[zone_idx]
        for f_axis_range_left_idx in range(len(zone_freqs) - 1):
            freq_left = zone_freqs[f_axis_range_left_idx]
            freq_right = zone_freqs[f_axis_range_left_idx + 1]
            imp_left_up = zone_upper_imps[f_axis_range_left_idx]
            imp_right_up = zone_upper_imps[f_axis_range_left_idx + 1]
            imp_left_down = zone_lower_imps[f_axis_range_left_idx]
            imp_right_down = zone_lower_imps[f_axis_range_left_idx + 1]
            if freq_left <= freq <= freq_right: #range of frequencies
                distance_out = (log_freq - math.log10(freq_left)) / (math.log10(freq_right) - math.log10(freq_left)) #calculates the distance away from the left endpoint the point is in relative terms
                upper_bound_imp = 10 ** (math.log10(imp_left_up) - ((math.log10(imp_left_up) - math.log10(imp_right_up)) * distance_out)) # calculates upper limit for impedance at the frequency you're at
                lower_bound_imp = 10 ** (math.log10(imp_left_down) - ((math.log10(imp_left_down) - math.log10(imp_right_down)) * distance_out)) # calculates lower limit for impedance at the frequency you're at
                if lower_bound_imp <= total_imp <= upper_bound_imp:
                    z_unc = z_unc_of_zone[zone_idx] #uncertainty we assign for this zone in length of vector in percent
                    ang_unc = ang_unc_of_zone[zone_idx] #uncertianty in argument in degrees
    
    re_left_unc =  0#gives the uncertainty as the lower bound of the real impedance value
    re_right_unc = 0#gives the uncertainty as the upper bound of the real impedance value
    re_max_unc = 0#gives the maximum uncertainty for range of the real impedance value
    neg_im_bottom_unc = 0#gives the uncertainty as the lower bound of the negative imaginary impedance value, or the upper bound of imaginary impedance
    neg_im_top_unc = 0#gives the uncertainty as the upper bound of the negative imaginary impedance value, or the lower bound of imaginary impedance
    neg_im_max_unc = 0#gives the maximum uncertainty for range of possible imaginary impedance values
    if ang_unc == 180:
        re_left_unc =  10 * total_imp
        re_right_unc = 10 * total_imp
        re_max_unc = 10 * total_imp
        neg_im_bottom_unc = 10 * total_imp
        neg_im_top_unc = 10 * total_imp
        neg_im_max_unc = 10 * total_imp
    else:
        ang_unc_rad = ang_unc * 2 * math.pi / 360
        re_corner_locations = []
        re_corner_locations.append((re*math.cos(ang_unc_rad) - neg_im * math.sin(ang_unc_rad)) * (100 + z_unc) / 100) #calculates the rotation with the length change
        re_corner_locations.append((re*math.cos(ang_unc_rad) - neg_im * math.sin(ang_unc_rad)) * (100 - z_unc) / 100)
        re_corner_locations.append((re*math.cos(-ang_unc_rad) - neg_im * math.sin(-ang_unc_rad)) * (100 + z_unc) / 100)
        re_corner_locations.append((re*math.cos(-ang_unc_rad) - neg_im * math.sin(-ang_unc_rad)) * (100 - z_unc) / 100)
        neg_im_corner_locations = []
        neg_im_corner_locations.append((re*math.sin(ang_unc_rad) + neg_im * math.cos(ang_unc_rad)) * (100 + z_unc) / 100)
        neg_im_corner_locations.append((re*math.sin(ang_unc_rad) + neg_im * math.cos(ang_unc_rad)) * (100 - z_unc) / 100)
        neg_im_corner_locations.append((re*math.sin(-ang_unc_rad) + neg_im * math.cos(-ang_unc_rad)) * (100 + z_unc) / 100)
        neg_im_corner_locations.append((re*math.sin(-ang_unc_rad) + neg_im * math.cos(-ang_unc_rad)) * (100 - z_unc) / 100)
        re_left_unc =  abs(min(re_corner_locations) - re)
        re_right_unc = abs(max(re_corner_locations) - re)
        re_max_unc = max(re_left_unc, re_right_unc)
        neg_im_bottom_unc = abs(min(neg_im_corner_locations) - neg_im)
        neg_im_top_unc = abs(max(neg_im_corner_locations) - neg_im)
        neg_im_max_unc = max(neg_im_bottom_unc, neg_im_top_unc)
    return [re_left_unc, re_right_unc, re_max_unc, neg_im_bottom_unc, neg_im_top_unc, neg_im_max_unc]
#____________________________error determiner_______

#make another function to convert frequency id to frequency


resistance_values_right_scan = [] # creates a list of resistance values found by scan for minimum from right
resistance_values_right_scan_error_low = []
resistance_values_right_scan_error_high = []
resistance_values_valley_scan = [] # creates a list of resistance values found by scan from bottom
date_and_time = [] # creates a list of the corresponding dates and times of data collection
value_counters_right_scan = [] # creates a list of which data point from the right scan was counted as the valley, indexed from 0, low to high freq, low to high R (generally)
value_counters_valley_scan = [] # creates a list of which data point from the bottom scan was counted as the valley, indexed from 0, low to high freq, low to high R (generally)
problematic_sets = [] # creates a list of data sets that are flagged as problematic

nyquist_indices = which_to_nyquist_plot #creates a list of which data sets were nyquist plotted
for ramp_idx in ramp_data_sets:
    if ramp_idx in nyquist_indices:
        nyquist_indices.remove(ramp_idx)
nyquist_plots = [] #creates a list of the nyquist plot data values, of the format [[[z re], [- z im]], [[z re], [- z im]], [[z re], [- z im]], ...], and z re is all of the real impedance values in a data set and so on
nyquist_file_name_splits = [] # creates a list of file name splits for naming graphs in nyquist plotting

for range_lst in plot_ranges_overplot:
    for idx in range_lst:
        if not (idx in nyquist_indices):
            nyquist_indices.append(idx)
            nyquist_indices.sort()

adjusted_times_minutes = [] #creates a list of adjusted times, in minutes
adjusted_times_seconds = [] #creates a list of times in seconds

ramp_dict = {} #TODO take this out once multiple rampins sets are accounted for
for i in range(1, (n_data_sets_total + 1)): # this is to scan through the files that you want to check. here, I scan through the second file to the 424th file
    file_name = f"C:/Users/Natha/Downloads/Coe College Research/EIS Data/2024-07-30 17H06M33S Nathan 25 Li2O 75 B2O3/2024-07-30 17H06M33S Nathan 25 Li2O 75 B2O3 - {i} of 2502.txt" # TODO creates a string specifying the file location
    data_file = open(file_name, "r") #opens the file
    if i in ramp_data_sets: # TODO account for multiple ramp data sets, just make another dictionary in front so we have a dictionary of dictionary of lists of pair lists
        line_counter = 0
        ramp_dict = {} #creates a dictionary of the form {temp1: [[z', -z''], [z', -z''],..., [z', -z'']], temp2: [[z', -z''], [z', -z''],..., [z', -z'']]...}
        for line in data_file: 
            line_split = line.split(",") #splits the values at the commas
            if line_counter >=4: #takes from the fourth line (indexed from zero), which was desired when writing this
                temp_ramp = float(line_split[ramp_data_sets_temp_column].strip()) #extracts temperature value
                z_real_ramp = float(line_split[ramp_data_sets_real_column].strip()) #extracts real impedance value
                neg_z_im_ramp = (-1 * float(line_split[ramp_data_sets_imaginary_column].strip())) #extracts negative imaginary impedance value
                re_neg_im_pair = [z_real_ramp, neg_z_im_ramp]
                if temp_ramp in ramp_dict.keys():
                    ramp_dict[temp_ramp].append(re_neg_im_pair) #adds the values to each ramp data set
                else:
                    ramp_dict[temp_ramp] = []
                    ramp_dict[temp_ramp].append(re_neg_im_pair)
            line_counter += 1
    else:
        line_counter = 0 # counts the line number
        z_real = [] # creates a list of real impedance values, in the order they came in, which is from the left
        neg_z_imaginary = [] # creates a list of imaginary impedance values, in the order they came in, which is from the left
        freq_idx = [] #creates a list of frequency indices, indexed from 0
        freq = [] #creates a list of frequencies, beginning from 10^7 or whatever the first one is
        unc_list = [] #list of list of uncertainty values as defined in the function above of the form [[re_left_unc, re_right_unc, re_max_unc, neg_im_bottom_unc, neg_im_top_unc, neg_im_max_unc], [re_left_unc, re_right_unc, re_max_unc, neg_im_bottom_unc, neg_im_top_unc, neg_im_max_unc], [re_left_unc, re_right_unc, re_max_unc, neg_im_bottom_unc, neg_im_top_unc, neg_im_max_unc]...]
        unc_bounds = [] #list of upper and lower bounds for each data point's two values of re and neg_im impedance of the form [[re_left_lower, re_right_upper, neg_im_bottom_lower, neg_im_top_upper], []...]
        for line in data_file: #extracts the desired values
            line_split = line.split(",") #splits the values at the commas
            if line_counter >=4 and line_counter < 4 + n_data_points_per_set: #takes the 4th through last, which was desired when I wrote this
                #print(line_split[3]) #FOR TEST
                temporary_neg_z_imaginary = -1 * float((line_split[imaginary_column]).strip())
                temporary_z_real = float((line_split[real_column]).strip())
                temporary_freq = float((line_split[freq_column]).strip())
                neg_z_imaginary.append(temporary_neg_z_imaginary)
                z_real.append(temporary_z_real)
                freq_idx.append(line_counter - 4)
                freq.append(temporary_freq) # with freq_column usually being the first column, with frequency values
                temprary_unc_list_entry_list_of_uncertainties_for_this_point = uncertainty_re_im(temporary_z_real, temporary_neg_z_imaginary, temporary_freq)
                unc_list.append(temprary_unc_list_entry_list_of_uncertainties_for_this_point)
                unc_bounds.append([temporary_z_real - temprary_unc_list_entry_list_of_uncertainties_for_this_point[0], temporary_z_real + temprary_unc_list_entry_list_of_uncertainties_for_this_point[1], temporary_neg_z_imaginary - temprary_unc_list_entry_list_of_uncertainties_for_this_point[3], temporary_neg_z_imaginary + temprary_unc_list_entry_list_of_uncertainties_for_this_point[4]])
            elif line_counter == 0: #takes the date and time
                date_and_time_split_temp = line_split[1:]
                date_and_time.append(line_split[1:]) #adds the date and time
            line_counter += 1
        data_file.close #closes the file VERY IMPORTANT

        #now, determine uncertainties

        #z_real.reverse() #reverses the data, so we can scan from the right
        #neg_z_imaginary.reverse()  #reverses the data, so we can scan from the right
        #print(z_real)  #FOR TEST
        #print(neg_z_imaginary)  #FOR TEST
        hit_local_minimum = False #status of whether you've reached the local minimum
        value_counter = -1 #sets a counter for which value you're on. we are doing reverse indexing
        flagged_value = intended_flagged_value - n_data_points_per_set # calculates the reverse indexing value
        forbidden_value_right = intended_forbidden_value_right - n_data_points_per_set # calculates the reverse indexing value
        forbidden_value_left = intended_forbidden_value_left - n_data_points_per_set # calculates the reverse indexing value
        number_above = 0 #creates a variable for the number of points that have been increasing, as you scan from the right


        #scanning for local minimum from the right
        while hit_local_minimum == False and value_counter > -n_data_points_per_set + 1: #we have n_data_points_per_set data points, and we want to compare pair, until we hit the latest possible value. however, we have a big problem if we hit the latest possible value
            if neg_z_imaginary[value_counter] > neg_z_imaginary[value_counter - 1]:
                number_above = 0
            elif neg_z_imaginary[value_counter] < neg_z_imaginary[value_counter - 1]: # starts increasing again, when scanning leftwards
                if (value_counter <= forbidden_value_right) and (value_counter >= forbidden_value_left): # if to the left of the maximum allowed value, then you are allowed
                    #resistance_values_right_scan.append(z_real[value_counter])
                    #value_counters_right_scan.append(value_counter + n_data_points_per_set)
                    #if value_counter >= flagged_value: #if to the right of the problematic index minimum
                    #    problematic_sets.append(i) #appends problematic index
                    #hit_local_minimum = True
                    number_above = number_above + 1
            elif neg_z_imaginary[value_counter] == neg_z_imaginary[value_counter - 1]: # catches edge case of a tie
                if (value_counter <= forbidden_value_right) and (value_counter >= forbidden_value_left): # if to the left of the maximum allowed value, then you are allowed
                    #resistance_values_right_scan.append(((z_real[value_counter] + z_real[value_counter - 1]) / 2))
                    #value_counters_right_scan.append(value_counter + n_data_points_per_set)
                    #if value_counter >= flagged_value: #if to the right of the problematic index minimum
                    #    problematic_sets.append(i) #appends problematic index
                    #hit_local_minimum = True
                    number_above = number_above + 1
            if number_above == min_rise_dist:
                resistance_values_right_scan.append(z_real[value_counter - 1 + min_rise_dist])
                value_counters_right_scan.append(value_counter - 1 + min_rise_dist + n_data_points_per_set)
                if value_counter >= flagged_value: #if to the right of the problematic index minimum
                    problematic_sets.append(i) #appends problematic index
                hit_local_minimum = True

                #uncertainty determination
                upper_bound_neg_im_minimum = unc_bounds[value_counter - 1 + min_rise_dist][3]
                minimum_idx = value_counter - 1 + min_rise_dist + n_data_points_per_set
                candidate_low_idxs_valley = [minimum_idx] #creates a list of candidate true minimum values
                hit_edge_pos = False #hit the last possible minimum value in to the right?
                hit_edge_neg = False #hit the last possible minimum value in to the right?
                steps_right = 0
                steps_left = 0
                while hit_edge_pos == False:
                    steps_right += 1
                    if (minimum_idx + steps_right < n_data_points_per_set) and (unc_bounds[minimum_idx + steps_right][2] < upper_bound_neg_im_minimum):
                        candidate_low_idxs_valley.append((minimum_idx + steps_right))
                    else:
                        hit_edge_pos = True
                while hit_edge_neg == False:
                    steps_left += 1
                    if (minimum_idx - steps_left >= 0) and (unc_bounds[minimum_idx - steps_left][2] < upper_bound_neg_im_minimum):
                        candidate_low_idxs_valley.append((minimum_idx - steps_left))
                    else:
                        hit_edge_neg = True
                for cand_min_idx in candidate_low_idxs_valley: #for each candidate true minimum in discrete case, the continuous fit might reach a minimum in between one such point and an adjacent point
                    if (cand_min_idx - 1 >= 0) and (cand_min_idx + 1 < len(z_real)) and (unc_bounds[cand_min_idx - 1][2] > unc_bounds[cand_min_idx + 1][3]): #if the low of the point with one less index is bigger than the high of the other adjacent point point, definitely the second lowest neg im point is the latter
                        if (cand_min_idx + 1) not in candidate_low_idxs_valley:
                            candidate_low_idxs_valley.append((cand_min_idx + 1))
                    elif (cand_min_idx - 1 >= 0) and (cand_min_idx + 1 < len(z_real)) and (unc_bounds[cand_min_idx - 1][3] < unc_bounds[cand_min_idx + 1][2]): #if the high of the point with one less index is smaller than the low of the other adjacent point point, definitely the second lowest neg im point is the first point
                        if (cand_min_idx - 1) not in candidate_low_idxs_valley:
                            candidate_low_idxs_valley.append((cand_min_idx - 1))
                    else:
                        if (cand_min_idx + 1 < len(z_real)) and ((cand_min_idx + 1) not in candidate_low_idxs_valley):
                            candidate_low_idxs_valley.append((cand_min_idx + 1))
                        if (cand_min_idx - 1 >= 0) and ((cand_min_idx - 1) not in candidate_low_idxs_valley):
                            candidate_low_idxs_valley.append((cand_min_idx - 1))
                min_cand_idx = min(candidate_low_idxs_valley)
                max_cand_idx = max(candidate_low_idxs_valley)
                dipping_point_low_bound = min([unc_bounds[min_cand_idx][0], unc_bounds[max_cand_idx][0]])
                dipping_point_high_bound = max([unc_bounds[min_cand_idx][1], unc_bounds[max_cand_idx][1]])
                resistance_values_right_scan_error_low.append(abs(z_real[value_counter - 1 + min_rise_dist] - dipping_point_low_bound)) #z_real[value_counter - 1 + min_rise_dist]
                resistance_values_right_scan_error_high.append(abs(z_real[value_counter - 1 + min_rise_dist] - dipping_point_high_bound)) #z_real[value_counter - 1 + min_rise_dist]

            value_counter = value_counter - 1 #reverse indexing step by one
        
        if hit_local_minimum == False:
            resistance_values_right_scan.append(-1)
            value_counters_right_scan.append(value_counter + min_rise_dist + n_data_points_per_set + intended_forbidden_value_left)
        
        #merge z_real, neg_z_imaginary, and position value
        merged_list = [] # creates a list of the form [[-z'', z', index], [-z'', z', index], ...] out of all three data sets
        for j in range(intended_forbidden_value_left, intended_forbidden_value_right + 1):
            merged_list.append([neg_z_imaginary[j], z_real[j], j])
        
        # print(merged_list) FOR TEST
        merged_sorted_list = sorted(merged_list)
        # print(merged_sorted_list) FOR TEST
        candidate_lows = [] #creates candidate list of minimum values with index of the point in the data set
        cut_candidate_lows = [] #creates sharper list of candidate lows

        #scan for second or only valley
        found_proper_valley = False
        k1 = 0
        while found_proper_valley == False and k1 <= intended_forbidden_value_right - intended_forbidden_value_left: #ends once you get to the limits of your index, maximum value, adds one by one to candidate lows and checks separation
            candidate_lows.append(merged_sorted_list[k1][2]) #adds the index of the point with lowest z''
            candidate_lows.sort() #sorts this list to put lowest index first
            for k2 in range(len(candidate_lows)-1): # checks to see if there is ever a gap of length low_separation between multiple groups of low values
                if abs(candidate_lows[k2] - candidate_lows[k2+1]) >= low_separation: # if any step between two values is greater than or equal to the separation
                    k2p1 = k2 + 1
                    cut_candidate_lows = candidate_lows[k2p1:] # slice off the first set
                    found_proper_valley = True
            k1 += 1

        if found_proper_valley == False:
            cut_candidate_lows = candidate_lows
        
        cut_merged_list = [] #all candidate data points, containing z'', z', and idx
        for k3 in cut_candidate_lows:
            cut_merged_list.append(merged_list[k3 - intended_forbidden_value_left]) #add back in the full objects
        valley_minimum = min(cut_merged_list) #take the smallest
        resistance_values_valley_scan.append(valley_minimum[1]) #append the resistance
        value_counters_valley_scan.append(valley_minimum[2]) #append the index

        if i == 1 + n_data_sets_exclude_at_start:
            adjusted_times_minutes.append(0) #if it's the first time, take the time to be 0
            current_time_seconds = 0
            adjusted_times_seconds.append(current_time_seconds)
        else:
            if i > 1 + n_data_sets_exclude_at_start:
                if date_and_time[i - (1 + n_data_sets_exclude_at_start) - 1][0] == date_and_time_split_temp[0]: #if the days are the same
                    time_lst_1 = date_and_time_split_temp[1].strip().split(":") #split into a list of hour and minutes
                    time_lst_2 = date_and_time[i - (1 + n_data_sets_exclude_at_start) - 1][1].strip().split(":")
                    time_1 = int(time_lst_1[0]) * 60 + int(time_lst_1[1]) #convert to minutes
                    time_2 = int(time_lst_2[0]) * 60 + int(time_lst_2[1])
                    previous_time = adjusted_times_minutes[i - (1 + n_data_sets_exclude_at_start) - 1] #takes the time of the index before TODO check
                    adjusted_times_minutes.append(previous_time + (time_1 - time_2)) #adds new time step to the previous time to get the current time
                    current_time_seconds = (previous_time + (time_1 - time_2)) * 60
                    adjusted_times_seconds.append(current_time_seconds)
                else: #if the days are different
                    time_lst_1 = date_and_time_split_temp[1].strip().split(":") #same as above
                    time_lst_2 = date_and_time[i - (1 + n_data_sets_exclude_at_start) - 1][1].strip().split(":")
                    time_1 = (int(time_lst_1[0]) + 24) * 60 + int(time_lst_1[1]) #same as above, but assumes that the later time is on the next day, so adds 24 hours to the hour time
                    time_2 = int(time_lst_2[0]) * 60 + int(time_lst_2[1])
                    previous_time = adjusted_times_minutes[i - (1 + n_data_sets_exclude_at_start) - 1] #takes the time of the index before TODO check
                    adjusted_times_minutes.append(previous_time + (time_1 - time_2)) #same as above
                    current_time_seconds = (previous_time + (time_1 - time_2)) * 60
                    adjusted_times_seconds.append(current_time_seconds)

        for time_range_lst in range_to_lin_fit_cond: #plots the nyquist plots of every plot within the range that one is taking the linear fit over
            if current_time_seconds >= time_range_lst[0] and current_time_seconds <= time_range_lst[1]: #checks if the time value is within one of the fitting ranges
                if not (i in nyquist_indices): #avoids duplication of nyquist plots
                    nyquist_indices.append(i)
                    nyquist_indices.sort()

        
        if i in nyquist_indices: #saves the data and a piece of the file name for any wanted nyquist plot in a form that's easy to read later
            nyquist_data_single = [] # of the form [[z re], [- z im]]
            nyquist_data_single.append(z_real)
            nyquist_data_single.append(neg_z_imaginary)
            file_name_split_right = file_name.split("/")[-1]
            nyquist_file_name_splits.append(file_name_split_right)
            nyquist_plots.append(nyquist_data_single) # TODO put the time in here as well



print("\n\n\nresistance values right scan:")
print(resistance_values_right_scan)
print(len(resistance_values_right_scan))
print("value counters right scan:")
print(value_counters_right_scan)
print("resistance values valley scan:")
print(resistance_values_valley_scan)
print(len(resistance_values_valley_scan))
print("value counters valley scan:")
print(value_counters_valley_scan)
print("dates and times:")
print(date_and_time)
print("problematic sets")
print(problematic_sets)

agreed_resistances = [] #creates a list of agreed-upon resistances
agreed_data_set_indices = [] # creates a list of the indices of the data sets that had agreed indices
agreed_resistance_errors_low = []
agreed_resistance_errors_high = []
for i2 in range(0, n_data_sets_total - n_data_sets_exclude_at_start):
    if resistance_values_right_scan[i2] == resistance_values_valley_scan[i2]:
        agreed_resistances.append(resistance_values_right_scan[i2])
        agreed_data_set_indices.append(i2 + 1) # indexes from one to match the file names
        agreed_resistance_errors_low.append(resistance_values_right_scan_error_low[i2])
        agreed_resistance_errors_high.append(resistance_values_right_scan_error_high[i2])

print("agreed resistances")
print(agreed_resistances)
print("agreed resistance rrors low")
print(agreed_resistance_errors_low)
print("agreed resistance errors high")
print(agreed_resistance_errors_high)
print("agreed data set indices")
print(agreed_data_set_indices)
print("overall times minutes")
print(adjusted_times_minutes)


print("adjusted times seconds")
print(adjusted_times_seconds)


agreed_times_seconds = [] #creates a list of agreed upon times in seconds
for i6 in agreed_data_set_indices:
    agreed_times_seconds.append(adjusted_times_seconds[i6 - 1])

print("agreed times seconds")
print(agreed_times_seconds)

agreed_resistivities = [] # creates a list of agreed upon resistivities, in ohm meters
agreed_conductivities = [] # creates a list of agreed upon conductivities, in ohm^-1 meters^-1
for resistance1 in agreed_resistances: #calculates resistivity and conductivity from resistance
    resistivity1 = (resistance1 * sample_electrode_area_mm2) / (sample_thickness_mm * 1000)
    agreed_resistivities.append(resistivity1)
    agreed_conductivities.append(1 / resistivity1)

print("agreed resistivities")
print(agreed_resistivities)
print("agreed conductivities")
print(agreed_conductivities)

data_to_fit_to = [] #creates a list of lists of lists of data that matlab will fit the overall resistance and conductivity values vs time to, of the form [[[seconds], [cond], [resistance]], [[seconds], [cond], [resistance]]...], with every [[seconds], [cond], [resistance]] being for a separate range of fitting
for range_lst in range_to_lin_fit_cond:
    fit_seconds_lst = []
    fit_conductivities_lst = []
    fit_resistances_lst = []
    for i6 in range(len(agreed_times_seconds)):
        if agreed_times_seconds[i6] >= range_lst[0] and agreed_times_seconds[i6] <= range_lst[1]:
            fit_seconds_lst.append(agreed_times_seconds[i6])
            fit_conductivities_lst.append(agreed_conductivities[i6])
            fit_resistances_lst.append(agreed_resistances[i6])
    data_to_fit_to.append([fit_seconds_lst, fit_conductivities_lst, fit_resistances_lst])

rejected_seconds = [] # creates a list of the seconds timings of the rejected data sets
for adj_sec in adjusted_times_seconds:
    if not (adj_sec in  agreed_times_seconds):
        rejected_seconds.append(adj_sec)
    
print("rejected seconds:")
print(rejected_seconds)
print("electrode dimensions:")
print("thickness (mm):")
print(sample_thickness_mm)
print("area (mm^2):")
print(sample_electrode_area_mm2)
if len(agreed_conductivities) >= 1:
    print("final conductivity value:")
    print(agreed_conductivities[-1])


print("MATLAB plotting code:\n\n\n") #prints formatted matlab code to plot this data
date_of_measurement = "_".join(((file_name.split("/")[-1]).split(" ")[0]).split("-"))
file_name_split_right = file_name.split("/")[-1]
date_of_measurement_spaces = " ".join(date_of_measurement.split("_"))

for temp_of_ramp in ramp_dict.keys():
    z_real_ramp_lst = []
    neg_z_im_ramp_lst = []
    for pair_lst in ramp_dict[temp_of_ramp]:
        z_real_ramp_lst.append(pair_lst[0])
        neg_z_im_ramp_lst.append(pair_lst[1])
    temp_of_ramp_name = str(temp_of_ramp).split('.')[0]
    print(f"tempRamp_z_real_{temp_of_ramp_name} = {z_real_ramp_lst}")
    print(f"tempRamp_neg_z_im_{temp_of_ramp_name} = {neg_z_im_ramp_lst}")
    print(f"plot(tempRamp_z_real_{temp_of_ramp_name}, tempRamp_neg_z_im_{temp_of_ramp_name}, 'k.')")
    print(f"title('Temperature of Ramp: {temp_of_ramp_name}')")


print("\n\n\n")
print(f"agreedTimesSeconds_{date_of_measurement} = {agreed_times_seconds}")
print(f"agreedConductivities_{date_of_measurement} = {agreed_conductivities}")
print(f"plot(agreedTimesSeconds_{date_of_measurement}, agreedConductivities_{date_of_measurement}, 'b.')")
print("hold on")

for i7 in range(len(data_to_fit_to)):
    print(f"agreed_seconds_cond_fit_{i7} = {data_to_fit_to[i7][0]}")
    print(f"agreed_conductivities_cond_fit_{i7} = {data_to_fit_to[i7][1]}")
    print(f"cond_fit_{i7} = fit(agreed_seconds_cond_fit_{i7}', agreed_conductivities_cond_fit_{i7}', 'poly1')")
    print(f"slope_cond_fit_{i7} = cond_fit_{i7}.p1 % ohms^{{-1}} meters^{{-1}} seconds^{{-1}}")
    print(f"predicted_conductivities_cond_fit_{i7} = cond_fit_{i7}.p2 + cond_fit_{i7}.p1 .* agreed_seconds_cond_fit_{i7}")
    print(f"plot(agreed_seconds_cond_fit_{i7}, predicted_conductivities_cond_fit_{i7})")

print("set(gca, 'fontsize', 20)")
print(f"title({{'{sample_name} heated at {heat_temperature}C ({heat_number} heat)', '{date_of_measurement_spaces}', '{file_name_split_right}'}})") # which heat
print("xlabel('Seconds of heating (Seconds)')")
print("ylabel('Equivalent Circuit Conductivity (Ohms^{-1} meters^{-1})')")
# print("axis([0 58000 0 .0012])") #optional to fix the plotting window
print("hold off\n\n\n")



print(f"agreedResistances_{date_of_measurement} = {agreed_resistances}")
print(f"resistanceErrorsNeg_{date_of_measurement} = {agreed_resistance_errors_low}")
print(f"resistanceErrorsPos_{date_of_measurement} = {agreed_resistance_errors_high}")
print(f"time_errors_{date_of_measurement} = 30 .* ones([1, {len(agreed_times_seconds)}])")
print(f"errorbar(agreedTimesSeconds_{date_of_measurement}, agreedResistances_{date_of_measurement}, resistanceErrorsNeg_{date_of_measurement}, resistanceErrorsPos_{date_of_measurement}, time_errors_{date_of_measurement}, time_errors_{date_of_measurement}, 'b.')")
print("hold on")

for i7 in range(len(data_to_fit_to)):
    print(f"agreed_seconds_res_fit_{i7} = {data_to_fit_to[i7][0]}")
    print(f"agreed_resistances_res_fit_{i7} = {data_to_fit_to[i7][2]}")
    print(f"res_fit_{i7} = fit(agreed_seconds_res_fit_{i7}', agreed_resistances_res_fit_{i7}', 'poly1')")
    print(f"slope_res_fit_{i7} = res_fit_{i7}.p1 % ohms^{{-1}} meters^{{-1}} seconds^{{-1}}")
    print(f"predicted_resistances_res_fit_{i7} = res_fit_{i7}.p2 + res_fit_{i7}.p1 .* agreed_seconds_res_fit_{i7}")
    print(f"plot(agreed_seconds_res_fit_{i7}, predicted_resistances_res_fit_{i7})")
    print(f"crystal_growth_rate = slope_res_fit_{i7} * {agreed_conductivities[-1]} * {sample_electrode_area_mm2} / 2000000") #m/s calculates growth rate of the crystal

print("set(gca, 'fontsize', 20)")
print(f"title({{'{sample_name} heated at {heat_temperature}C ({heat_number} heat)', '{date_of_measurement_spaces}', '{file_name_split_right}'}})") # which heat
print("xlabel('Seconds of heating (Seconds)')")
print("ylabel('Equivalent Circuit Resistance Ohms')")
#print("axis([0 65000 0 .00035])") #optional to fix the plotting window
print("hold off\n\n\n")



for i in range(len(nyquist_indices)):
    print(f"z_real_new_{adjusted_times_seconds[nyquist_indices[i] - (1 + n_data_sets_exclude_at_start)]} = ", nyquist_plots[i][0])
    print(f"neg_z_imaginary_{adjusted_times_seconds[nyquist_indices[i] - (1 + n_data_sets_exclude_at_start)]} = ", nyquist_plots[i][1])
    print(f"plot(z_real_new_{adjusted_times_seconds[nyquist_indices[i] - (1 + n_data_sets_exclude_at_start)]}, neg_z_imaginary_{adjusted_times_seconds[nyquist_indices[i] - (1 + n_data_sets_exclude_at_start)]}, 'k.')")
    #print("axis([0 250 0 250])") #optional line of code to fix the graph window
    print(f"title({{'Time in Seconds: {adjusted_times_seconds[nyquist_indices[i] - (1 + n_data_sets_exclude_at_start)]}', '{nyquist_file_name_splits[i]}'}})")


marker_shapes = ['o', '^', 'pentagram', '>', 'square', 'v', 'hexagram', '<', 'diamond']
for range_lst in plot_ranges_overplot:
    color_list = list(range(0, len(range_lst), 1)) #remember to divide by (len(range_lst) - 1) to normalize
    for i7 in range(len(range_lst)):
        if i7 == 0:
            print(f"\n\nplot(z_real_new_{adjusted_times_seconds[range_lst[i7] - (1 + n_data_sets_exclude_at_start)]}, neg_z_imaginary_{adjusted_times_seconds[range_lst[i7] - (1 + n_data_sets_exclude_at_start)]}, '{marker_shapes[(i7 % len(marker_shapes))]}', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [{color_list[i7] / (len(range_lst) - 1)}, 0, {color_list[-1 - i7] / (len(range_lst) - 1)}])")
            print("hold on")
        else:
            print(f"plot(z_real_new_{adjusted_times_seconds[range_lst[i7] - (1 + n_data_sets_exclude_at_start)]}, neg_z_imaginary_{adjusted_times_seconds[range_lst[i7] - (1 + n_data_sets_exclude_at_start)]}, '{marker_shapes[(i7 % len(marker_shapes))]}', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [{color_list[i7] / (len(range_lst) - 1)}, 0, {color_list[-1 - i7] / (len(range_lst) - 1)}])")
    print(f"legend(", end='')
    for i8 in range(len(range_lst)):
        if i8 != len(range_lst) - 1:
            print(f"'{adjusted_times_seconds[range_lst[i8] - (1 + n_data_sets_exclude_at_start)]} sec', ", end='')
        else:
            print(f"'{adjusted_times_seconds[range_lst[i8] - (1 + n_data_sets_exclude_at_start)]} sec')")
    print("hold off\n\n")


print("\n\n\n")
