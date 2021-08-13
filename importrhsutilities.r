load_file = function(filename) {
    # Start timing
    start = Sys.time()
    
    # Open file
    fid = file(filename, "rb")
    filesize = file.info(filename)$size
    
    # Read file header
    header = read_header(fid)
    
    # Output a summary of recorded data
    cat(paste("Found ", header$num_amplifier_channels, " amplifier channel", plural(header$num_amplifier_channels), "\n", sep = ""))
    cat(paste("Found ", header$num_board_adc_channels, " board ADC channel", plural(header$num_board_adc_channels), "\n", sep = ""))
    cat(paste("Found ", header$num_board_dac_channels, " board DAC channel", plural(header$num_board_dac_channels), "\n", sep = ""))
    cat(paste("Found ", header$num_board_dig_in_channels, " board digital input channel", plural(header$num_board_dig_in_channels), "\n", sep = ""))
    cat(paste("Found ", header$num_board_dig_out_channels, " board digital output channel", plural(header$num_board_dig_out_channels), "\n", sep = ""))
    
    # Determine how many samples the data file contains
    bytes_per_block = get_bytes_per_data_block(header)
    
    # Calculate how many data blocks are present
    data_present = 0
    bytes_remaining = filesize - seek(fid)
    if (bytes_remaining > 0) {
        data_present = 1
    }
    
    if ((bytes_remaining %% bytes_per_block) != 0) {
        stop("Something is wrong with file size: should have a whole number of data blocks")
    }
    
    num_data_blocks = bytes_remaining / bytes_per_block
    
    # Calculate how many samples of each signal type are present
    num_amplifier_samples = 128 * num_data_blocks
    num_board_adc_samples = 128 * num_data_blocks
    num_board_dac_samples = 128 * num_data_blocks
    num_board_dig_in_samples = 128 * num_data_blocks
    num_board_dig_out_samples = 128 * num_data_blocks
    
    # Calculate how much time has been recorded
    record_time = num_amplifier_samples / header$sample_rate
    
    # Output a summary of contents of header file
    if (data_present > 0) {
        cat(paste("File contains ", record_time, " seconds of data. Amplifiers were sampled at ", header$sample_rate / 1000, " kS/s.\n", sep = ""))
    } else {
        cat(paste("Header file contains no data. Amplifiers were sampled at ", header$sample_rate / 1000, " kS/s.\n", sep = ""))
    }
    
    if (data_present > 0) {
        # Pre-allocate memory for data
        cat("Allocating memory for data...\n")
        data = list()
        data$t = matrix(0, 1, num_amplifier_samples)
        
        data$amplifier_data = allocate_zeros(header$num_amplifier_channels, num_amplifier_samples)
        
        if (header$dc_amplifier_data_saved > 0) {
            data$dc_amplifier_data = allocate_zeros(header$num_amplifier_channels, num_amplifier_samples)
        }
        
        data$stim_data_raw = allocate_zeros(header$num_amplifier_channels, num_amplifier_samples)
        data$stim_data = allocate_zeros(header$num_amplifier_channels, num_amplifier_samples)
        data$compliance_limit_data = allocate_zeros(header$num_amplifier_channels, num_amplifier_samples)
        data$charge_recovery_data = allocate_zeros(header$num_amplifier_channels, num_amplifier_samples)
        data$amp_settle_data = allocate_zeros(header$num_amplifier_channels, num_amplifier_samples)
        data$stim_polarity = allocate_zeros(header$num_amplifier_channels, num_amplifier_samples)
        
        data$board_adc_data = allocate_zeros(header$num_board_adc_channels, num_board_adc_samples)
        data$board_dac_data = allocate_zeros(header$num_board_dac_channels, num_board_dac_samples)
        
        # By default, this script interprets digital events (digital inputs and outputs) as Int16
        data$board_dig_in_data = allocate_zeros(header$num_board_dig_in_channels, num_board_dig_in_samples)
        data$board_dig_in_raw = allocate_zeros(1, num_board_dig_in_samples)
        data$board_dig_out_data = allocate_zeros(header$num_board_dig_out_channels, num_board_dig_out_samples)
        data$board_dig_out_raw = allocate_zeros(1, num_board_dig_out_samples)
        
        # Read sampled data from file
        cat("Reading data from file...\n")
        
        # Initialize indices used in looping
        indices = list("amplifier" = 1)
        indices$board_adc = 1
        indices$board_dac = 1
        indices$board_dig_in = 1
        indices$board_dig_out = 1
        
        print_increment = 10
        percent_done = print_increment
        for (i in 1:num_data_blocks) {
            data = read_one_data_block(data, header, indices, fid)
            
            # Increment indices
            indices$amplifier = 128 + indices$amplifier
            indices$board_adc = 128 + indices$board_adc
            indices$board_dac = 128 + indices$board_dac
            indices$board_dig_in = 128 + indices$board_dig_in
            indices$board_dig_out = 128 + indices$board_dig_out
            
            fraction_done = 100 * (1.0 * i / num_data_blocks)
            if (fraction_done >= percent_done) {
                cat(paste(percent_done, "% done...\n", sep = ""))
                percent_done = percent_done + print_increment
            }
        }
        
        # Make sure we have read exactly the right amount of data
        bytes_remaining = filesize - seek(fid)
        if (bytes_remaining != 0) {
            stop("Error: End of file not reached.")
        }
    } else {
        data = NA
    }
    
    if (data_present > 0) {
        cat("Parsing data...\n")
        
        # Extract digital input channels to separate variables
        if (header$num_board_dig_in_channels > 0) {
            for (i in 1:header$num_board_dig_in_channels) {
                mask = 2^(header$board_dig_in_channels[[i]]$native_order)
                data$board_dig_in_data[i, 1:ncol(data$board_dig_in_data)] = extract_digital(data$board_dig_in_raw[1, 1:ncol(data$board_dig_in_data)], mask)
            }
        }
        
        # Extract digital output channels to separate variables
        if (header$num_board_dig_out_channels > 0) {
            for (i in 1:header$num_board_dig_out_channels) {
                mask = 2^(header$board_dig_out_channels[[i]]$native_order)
                data$board_dig_out_data[i, 1:ncol(data$board_dig_out_data)] = extract_digital(data$board_dig_out_raw[1, 1:ncol(data$board_dig_out_data)], mask)
            }
        }
        
        # Extract stimulation data
        if (header$num_amplifier_channels > 0) {
            for (i in 1:header$num_amplifier_channels) {
                compliance_mask = 32768 # 2^15 bit
                data$compliance_limit_data[i, 1:ncol(data$compliance_limit_data)] = extract_digital(data$stim_data_raw[i, 1:ncol(data$compliance_limit_data)], compliance_mask)
                
                charge_recovery_mask = 16384 # 2^14 bit
                data$charge_recovery_data[i, 1:ncol(data$charge_recovery_data)] = extract_digital(data$stim_data_raw[i, 1:ncol(data$charge_recovery_data)], charge_recovery_mask)
                
                amp_settle_mask = 8192 # 2^13 bit
                data$amp_settle_data[i, 1:ncol(data$amp_settle_data)] = extract_digital(data$stim_data_raw[i, 1:ncol(data$amp_settle_data)], amp_settle_mask)
                
                stim_polarity_mask = 256 # 2^8 bit
                data$stim_polarity[i, 1:ncol(data$stim_polarity)] = 1 - (2 * extract_digital(data$stim_data_raw[i, 1:ncol(data$stim_polarity)], stim_polarity_mask))
                
                curr_amp = bitwAnd(data$stim_data_raw[i, 1:ncol(data$stim_data)], rep(255, length(data$stim_data_raw[i, 1:ncol(data$stim_data)]))) # Magnitude is lower 8 bits
                
                data$stim_data[i, 1:ncol(data$stim_data)] = curr_amp * data$stim_polarity[i, 1:ncol(data$stim_polarity)] # multiply current amplitude by the correct sign
            }
        }
                
        # Scale voltage levels appropriately
        data$amplifier_data = 0.195 * (data$amplifier_data - 32768) # units = microvolts
        data$stim_data = header$stim_step_size * (data$stim_data / 1.0e-6)
        if (header$dc_amplifier_data_saved > 0) {
            data$dc_amplifier_data = -0.01923 * (data$dc_amplifier_data - 512) # units = volts
        }
        data$board_adc_data = 312.5e-6 * (data$board_adc_data - 32768) # units = volts
        data$board_dac_data = 312.5e-6 * (data$board_dac_data - 32768) # units = volts

        # Check for gaps in timestamps
        num_gaps = 0

        diff_vector = diff(c(data$t))
        for (i in 1:length(diff_vector)) {
            if (diff_vector[i] != 1) {
                num_gaps = num_gaps + 1
            }
        }

        if (num_gaps == 0) {
            cat("No missing timestamps in data.\n")
        } else {
            cat(paste("Warning: ", num_gaps, " gaps in timestamp data found. Time scale will not be uniform!\n", sep = ""))
        }

        # Scale time steps (units = seconds).
        data$t = data$t / header$sample_rate

        # If the software notch filter was selected during the recording, apply the same notch filter to amplifier data here
        if (header$frequency_parameters$notch_filter_frequency > 0 && header$version$major < 3) {
            cat("Applying notch filter...\n")

            print_increment = 10
            percent_done = print_increment
            for (i in 1:header$num_amplifier_channels) {
                data$amplifier_data[i, 1:ncol(data$amplifier_data)] = notch_filter(data$amplifier_data[i, 1:ncol(data$amplifier_data)], header$sample_rate, header$frequency_parameters$notch_filter_frequency, 10)
                fraction_done = 100 * (i / header$num_amplifier_channels)
                if (fraction_done >= percent_done) {
                    cat(paste(percent_done, "% done...\n", sep = ""))
                    percent_done = percent_done + print_increment
                }
            }
        }
    }
    
    # Move variables to result variable
    result = data_to_result(header, data, data_present)
    
    # End timing
    elapsed = Sys.time() - start
    cat(paste("Done! Elapsed time: ", elapsed, " seconds\n", sep = ""))
    
    return(list("result" = result, "data_present" = data_present))
}

# Define print_all_channel_names function
print_all_channel_names = function(result) {
    # Print all amplifier_channels
    print_names_in_group(result$amplifier_channels)
    
    # Print all dc_amplifier_channels
    if (result$dc_amplifier_data_saved > 0) {
        print_names_in_group(result$dc_amplifier_channels)
    }
    
    # Print all stim_channels
    print_names_in_group(result$stim_channels)
    
    # Print all amp_settle_channels
    print_names_in_group(result$amp_settle_channels)
    
    # Print all charge_recovery_channels
    print_names_in_group(result$charge_recovery_channels)
    
    # Print all compliance_limit_channels
    print_names_in_group(result$compliance_limit_channels)
    
    # Print all board_adc_channels
    print_names_in_group(result$board_adc_channels)
    
    # Print all board_dac_channels
    print_names_in_group(result$board_dac_channels)
    
    # Print all board_dig_in_channels
    print_names_in_group(result$board_dig_in_channels)
    
    # Print all board_dig_out_channels
    print_names_in_group(result$board_dig_out_channels)
}

# Define function print_names_in_group function
print_names_in_group = function(signal_group) {
    for (this_channel in signal_group) {
        cat(paste(this_channel$custom_channel_name, "\n", sep=""))
    }
}

extract_digital = function(raw_data, mask) {
    mask_vector = rep(mask, length(raw_data))
    intermediate = bitwAnd(raw_data, mask_vector)
    result = c()
    for (i in 1:length(raw_data)) {
        if (intermediate[i]) result[i] = 1 else result[i] = 0
    }
    
    return(result)
}

# Define allocate_zeros
allocate_zeros = function(num_channels, num_samples) {
    if (num_channels > 0) {
        return(matrix(0, num_channels, num_samples))
    } else {
        return(matrix(, nrow=0, ncol=0)) # return empty matrix
    }
}

# Define get_bytes_per_data_block
get_bytes_per_data_block = function(header) {
    # Calculates the number of bytes in each 128-sample datablock.
    N = 128 # n of amplifier samples
    
    # Each data block contains N amplifier samples
    bytes_per_block = N * 4 # timestamp data
    
    bytes_per_block = bytes_per_block + N * 2 * header$num_amplifier_channels
    
    # DC amplifier voltage (absent if flag was off)
    if (header$dc_amplifier_data_saved > 0) {
        bytes_per_block = bytes_per_block + N * 2 * header$num_amplifier_channels
    }
    
    # Stimulation data, one per enabled amplifier channel
    bytes_per_block = bytes_per_block + N * 2 * header$num_amplifier_channels
    
    # Board analog inputs are sampled at same rate as amplifiers
    bytes_per_block = bytes_per_block + N * 2 * header$num_board_adc_channels
    
    # Board analog outputs are sampled at same rate as amplifiers
    bytes_per_block = bytes_per_block + N * 2 * header$num_board_dac_channels
    
    # Board digital inputs are sampled at same rate as amplifiers
    if (header$num_board_dig_in_channels > 0) {
        bytes_per_block = bytes_per_block + N * 2
    }
    
    # Board digital outputs are sampled at same rate as amplifiers
    if (header$num_board_dig_out_channels > 0) {
        bytes_per_block = bytes_per_block + N * 2
    }
    
    return(bytes_per_block)
}

# Define read_header function
read_header = function(fid) {
    # Check 'magic number' at beginning of file to make sure this is an Intan Technologies RHS2000 data file.
    magic_number = readUInt32(fid)
    if (magic_number != 0xd69127ac) {
        stop("Unrecognized file type")
    }
    
    # Read version number
    version = list(
        major = readInt16(fid),
        minor = readInt16(fid)
    )
    header = list("version" = version)
    
    cat(paste("Reading Intan Technologies RHS2000 Data File, Version ", version$major, ".", version$minor, "\n", sep = ""))
    
    # Read information of sampling rate and amplifier frequency settings
    header$sample_rate = readFloat32(fid)
    
    frequency_parameters = list("dsp_enabled" = readInt16(fid))
    frequency_parameters$actual_dsp_cutoff_frequency = readFloat32(fid)
    frequency_parameters$actual_lower_bandwidth = readFloat32(fid)
    frequency_parameters$actual_lower_settle_bandwidth = readFloat32(fid)
    frequency_parameters$actual_upper_bandwidth = readFloat32(fid)
    frequency_parameters$desired_dsp_cutoff_frequency = readFloat32(fid)
    frequency_parameters$desired_lower_bandwidth = readFloat32(fid)
    frequency_parameters$desired_lower_settle_bandwidth = readFloat32(fid)
    frequency_parameters$desired_upper_bandwidth = readFloat32(fid)
    
    # This tells us if a software 50/60 Hz notch filter was enabled during the data acquisition
    notch_filter_mode = readInt16(fid)
    notch_filter_frequency = 0
    if (notch_filter_mode == 1) {
        notch_filter_frequency = 50
    } else if (notch_filter_mode == 2) {
        notch_filter_frequency = 60
    }
    frequency_parameters$notch_filter_frequency = notch_filter_frequency
    frequency_parameters$desired_impedance_test_frequency = readFloat32(fid)
    frequency_parameters$actual_impedance_test_frequency = readFloat32(fid)
    
    header$amp_settle_mode = readInt16(fid)
    header$charge_recovery_mode = readInt16(fid)
    
    header$stim_step_size = readFloat32(fid)
    header$recovery_current_limit = readFloat32(fid)
    header$recovery_target_voltage = readFloat32(fid)
    
    # Place notes in array of Strings
    header$notes = c(readQString(fid), readQString(fid), readQString(fid))
    
    header$dc_amplifier_data_saved = readInt16(fid)
    header$eval_board_mode = readInt16(fid)
    header$reference_channel = readQString(fid)
    
    # Place frequency-related information in data structure
    frequency_parameters$amplifier_sample_rate = header$sample_rate
    
    header$frequency_parameters = frequency_parameters
    
    header$spike_triggers = list()
    
    header$amplifier_channels = list()
    header$dc_amplifier_channels = list()
    header$stim_channels = list()
    header$amp_settle_channels = list()
    header$charge_recovery_channels = list()
    header$compliance_limit_channels = list()
    header$board_adc_channels = list()
    header$board_dac_channels = list()
    header$board_dig_in_channels = list()
    header$board_dig_out_channels = list()
    
    amplifier_index = 1
    board_adc_index = 1
    board_dac_index = 1
    board_dig_in_index = 1
    board_dig_out_index = 1
    
    # Read signal summary from data file header
    number_of_signal_groups = readInt16(fid)
    
    for (signal_group in 1:number_of_signal_groups) {
      
        signal_group_name = readQString(fid)
        signal_group_prefix = readQString(fid)
        signal_group_enabled = readInt16(fid)
        signal_group_num_channels = readInt16(fid)
        signal_group_num_amp_channels = readInt16(fid)
        
        if ((signal_group_num_channels > 0) && (signal_group_enabled > 0)) {
            for (signal_channel in 1:signal_group_num_channels) {
                
                new_channel = list("port_name" = signal_group_name)
                new_channel$port_prefix = signal_group_prefix
                new_channel$port_number = signal_group
                
                new_channel$native_channel_name = readQString(fid)
                new_channel$custom_channel_name = readQString(fid)
                new_channel$native_order = readInt16(fid)
                new_channel$custom_order = readInt16(fid)
                signal_type = readInt16(fid)
                channel_enabled = readInt16(fid)
                new_channel$chip_channel = readInt16(fid)
                command_stream = readInt16(fid)
                new_channel$board_stream = readInt16(fid)
                new_trigger_channel = list("voltage_trigger_mode" = readInt16(fid))
                new_trigger_channel$voltage_threshold = readInt16(fid)
                new_trigger_channel$digital_trigger_channel = readInt16(fid)
                new_trigger_channel$digital_edge_polarity = readInt16(fid)
                new_channel$electrode_impedance_magnitude = readFloat32(fid)
                new_channel$electrode_impedance_phase = readFloat32(fid)
                
                if (channel_enabled > 0) {
                    if (signal_type == 0) {
                        header$amplifier_channels = append(header$amplifier_channels, list(new_channel))
                        
                        # If dc amplifier data is being saved, dc_amplifier_channels
                        if (header$dc_amplifier_data_saved > 0) {
                            new_dc_channel = list()
                            new_dc_channel$native_channel_name = paste("DC_", new_channel$native_channel_name, sep="")
                            new_dc_channel$custom_channel_name = paste("DC_", new_channel$custom_channel_name, sep="")
                            new_dc_channel$native_order = new_channel$native_order
                            new_dc_channel$custom_order = new_channel$custom_order
                            new_dc_channel$board_stream = new_channel$board_stream
                            new_dc_channel$chip_channel = new_channel$chip_channel
                            new_dc_channel$port_name = new_channel$port_name
                            new_dc_channel$port_prefix = new_channel$port_prefix
                            new_dc_channel$port_number = new_channel$port_number
                            new_dc_channel$electrode_impedance_magnitude = new_channel$electrode_impedance_magnitude
                            new_dc_channel$electrode_impedance_phase = new_channel$electrode_impedance_phase
                            header$dc_amplifier_channels = append(header$dc_amplifier_channels, list(new_dc_channel))
                        }
                        
                        # stim_channels
                        new_stim_channel = list()
                        new_stim_channel$native_channel_name = paste("STIM_", new_channel$native_channel_name, sep="")
                        new_stim_channel$custom_channel_name = paste("STIM_", new_channel$custom_channel_name, sep="")
                        new_stim_channel$native_order = new_channel$native_order
                        new_stim_channel$custom_order = new_channel$custom_order
                        new_stim_channel$board_stream = new_channel$board_stream
                        new_stim_channel$chip_channel = new_channel$chip_channel
                        new_stim_channel$port_name = new_channel$port_name
                        new_stim_channel$port_prefix = new_channel$port_prefix
                        new_stim_channel$port_number = new_channel$port_number
                        new_stim_channel$electrode_impedance_magnitude = new_channel$electrode_impedance_magnitude
                        new_stim_channel$electrode_impedance_phase = new_channel$electrode_impedance_phase
                        header$stim_channels = append(header$stim_channels, list(new_stim_channel))
                        
                        # amp_settle_channels
                        new_amp_settle_channel = list()
                        new_amp_settle_channel$native_channel_name = paste("AMP_SETTLE_", new_channel$native_channel_name, sep="")
                        new_amp_settle_channel$custom_channel_name = paste("AMP_SETTLE_", new_channel$custom_channel_name, sep="")
                        new_amp_settle_channel$native_order = new_channel$native_order
                        new_amp_settle_channel$custom_order = new_channel$custom_order
                        new_amp_settle_channel$board_stream = new_channel$board_stream
                        new_amp_settle_channel$chip_channel = new_channel$chip_channel
                        new_amp_settle_channel$port_name = new_channel$port_name
                        new_amp_settle_channel$port_prefix = new_channel$port_prefix
                        new_amp_settle_channel$port_number = new_channel$port_number
                        new_amp_settle_channel$electrode_impedance_magnitude = new_channel$electrode_impedance_magnitude
                        new_amp_settle_channel$electrode_impedance_phase = new_channel$electrode_impedance_phase
                        header$amp_settle_channels = append(header$amp_settle_channels, list(new_amp_settle_channel))
                        
                        # charge_recovery_channels
                        new_charge_recovery_channel = list()
                        new_charge_recovery_channel$native_channel_name = paste("CHARGE_RECOVERY_", new_channel$native_channel_name, sep="")
                        new_charge_recovery_channel$custom_channel_name = paste("CHARGE_RECOVERY_", new_channel$custom_channel_name, sep="")
                        new_charge_recovery_channel$native_order = new_channel$native_order
                        new_charge_recovery_channel$custom_order = new_channel$custom_order
                        new_charge_recovery_channel$board_stream = new_channel$board_stream
                        new_charge_recovery_channel$chip_channel = new_channel$chip_channel
                        new_charge_recovery_channel$port_name = new_channel$port_name
                        new_charge_recovery_channel$port_prefix = new_channel$port_prefix
                        new_charge_recovery_channel$port_number = new_channel$port_number
                        new_charge_recovery_channel$electrode_impedance_magnitude = new_channel$electrode_impedance_magnitude
                        new_charge_recovery_channel$electrode_impedance_phase = new_channel$electrode_impedance_phase
                        header$charge_recovery_channels = append(header$charge_recovery_channels, list(new_charge_recovery_channel))
                        
                        # compliance_limit_channels
                        new_compliance_limit_channel = list()
                        new_compliance_limit_channel$native_channel_name = paste("COMPLIANCE_LIMIT_", new_channel$native_channel_name, sep="")
                        new_compliance_limit_channel$custom_channel_name = paste("COMPLIANCE_LIMIT_", new_channel$custom_channel_name, sep="")
                        new_compliance_limit_channel$native_order = new_channel$native_order
                        new_compliance_limit_channel$custom_order = new_channel$custom_order
                        new_compliance_limit_channel$board_stream = new_channel$board_stream
                        new_compliance_limit_channel$chip_channel = new_channel$chip_channel
                        new_compliance_limit_channel$port_name = new_channel$port_name
                        new_compliance_limit_channel$port_prefix = new_channel$port_prefix
                        new_compliance_limit_channel$port_number = new_channel$port_number
                        new_compliance_limit_channel$electrode_impedance_magnitude = new_channel$electrode_impedance_magnitude
                        new_compliance_limit_channel$electrode_impedance_phase = new_channel$electrode_impedance_phase
                        header$compliance_limit_channels = append(header$compliance_limit_channels, list(new_compliance_limit_channel))
                        
                        header$spike_triggers = append(header$spike_triggers, list(new_trigger_channel))
                        amplifier_index = amplifier_index + 1
                    }
                    else if (signal_type == 1) {
                        stop("Wrong signal type for the rhs format")
                    }
                    else if (signal_type == 2) {
                        stop("Wrong signal type for the rhs format")
                    }
                    else if (signal_type == 3) {
                        header$board_adc_channels = append(header$board_adc_channels, list(new_channel))
                        board_adc_index = board_adc_index + 1
                    }
                    else if (signal_type == 4) {
                        header$board_dac_channels = append(header$board_dac_channels, list(new_channel))
                        board_dac_index = board_dac_index + 1
                    }
                    else if (signal_type == 5) {
                        header$board_dig_in_channels = append(header$board_dig_in_channels, list(new_channel))
                        board_dig_in_index = board_dig_in_index + 1
                    }
                    else if (signal_type == 6) {
                        header$board_dig_out_channels = append(header$board_dig_out_channels, list(new_channel))
                        board_dig_out_index = board_dig_out_index + 1
                    }
                    else {
                        stop("Unknown channel type")
                    }
                }
            }
        }
    }
    
    # Summarize contents of data file
    header$num_amplifier_channels = amplifier_index - 1
    header$num_board_adc_channels = board_adc_index - 1
    header$num_board_dac_channels = board_dac_index - 1
    header$num_board_dig_in_channels = board_dig_in_index - 1
    header$num_board_dig_out_channels = board_dig_out_index - 1
    
    return(header)
}

# Define data_to_result function
data_to_result = function(header, data, data_present) {
    # Moves the header and data (if present) into a common object
    result = list()
    
    stim_parameters = list()
    stim_parameters$stim_step_size = header$stim_step_size
    stim_parameters$charge_recovery_current_limit = header$recovery_current_limit
    stim_parameters$charge_recovery_target_voltage = header$recovery_target_voltage
    stim_parameters$amp_settle_mode = header$amp_settle_mode
    stim_parameters$charge_recovery_mode = header$charge_recovery_mode
    result$stim_parameters = stim_parameters
    
    result$stim_channels = header$stim_channels
    result$spike_triggers = header$spike_triggers
    result$notes = header$notes
    result$frequency_parameters = header$frequency_parameters
    
    result$reference_channel = header$reference_channel
    result$amplifier_channels = header$amplifier_channels
    result$board_adc_channels = header$board_adc_channels
    result$board_dac_channels = header$board_dac_channels
    
    result$dc_amplifier_data_saved = header$dc_amplifier_data_saved
    if (header$dc_amplifier_data_saved > 0) {
        result$dc_amplifier_channels = header$dc_amplifier_channels
    }
    
    result$compliance_limit_channels = header$compliance_limit_channels
    result$charge_recovery_channels = header$charge_recovery_channels
    result$amp_settle_channels = header$amp_settle_channels
    result$board_dig_in_channels = header$board_dig_in_channels
    result$board_dig_out_channels = header$board_dig_out_channels
    
    if (data_present) {
        result$t = data$t
        result$stim_data = data$stim_data
        result$amplifier_data = data$amplifier_data
        result$board_adc_data = data$board_adc_data
        result$board_dac_data = data$board_dac_data
        if (header$dc_amplifier_data_saved > 0) {
            result$dc_amplifier_data = data$dc_amplifier_data
        }
        result$compliance_limit_data = data$compliance_limit_data
        result$charge_recovery_data = data$charge_recovery_data
        result$amp_settle_data = data$amp_settle_data
        result$board_dig_in_data = data$board_dig_in_data
        result$board_dig_out_data = data$board_dig_out_data
    }
    
    return(result)
}


# Define read_one_data_block function
read_one_data_block = function(data, header, indices, fid) {
    # Reads one 128 sample data block from fid into data, at the location indicated by indices
    data$t[indices$amplifier:(indices$amplifier + 128 - 1)] = readNInt32(fid, 128)
    
    if (header$num_amplifier_channels > 0) {
        data$amplifier_data[1:nrow(data$amplifier_data), indices$amplifier:(indices$amplifier + 128 - 1)] = matrix(readNUInt16(fid, nrow(data$amplifier_data) * 128), nrow=nrow(data$amplifier_data), byrow=TRUE)
    
        # Check if dc amplifier voltage was saved
        if (header$dc_amplifier_data_saved > 0) {
            data$dc_amplifier_data[1:nrow(data$dc_amplifier_data), indices$amplifier:(indices$amplifier + 128 - 1)] = matrix(readNUInt16(fid, nrow(data$dc_amplifier_data) * 128), nrow=nrow(data$dc_amplifier_data), byrow=TRUE)
        }
    
        # Get the stimulation data
        data$stim_data_raw[1:nrow(data$stim_data_raw), indices$amplifier:(indices$amplifier + 128 - 1)] = matrix(readNUInt16(fid, nrow(data$stim_data_raw) * 128), nrow=nrow(data$stim_data_raw), byrow=TRUE)
    }
    
    if (header$num_board_adc_channels > 0) {
        data$board_adc_data[1:nrow(data$board_adc_data), indices$board_adc:(indices$board_adc + 128 - 1)] = matrix(readNUInt16(fid, nrow(data$board_adc_data) * 128), nrow=nrow(data$board_adc_data), byrow=TRUE)
    }
    
    if (header$num_board_dac_channels > 0) {
        data$board_dac_data[1:nrow(data$board_dac_data), indices$board_dac:(indices$board_dac + 128 - 1)] = matrix(readNUInt16(fid, nrow(data$board_dac_data) * 128), nrow=nrow(data$board_dac_data), byrow=TRUE)
    }
    
    if (header$num_board_dig_in_channels > 0) {
        data$board_dig_in_raw[1, indices$board_dig_in:(indices$board_dig_in + 128 - 1)] = matrix(readNUInt16(fid, 128), nrow=nrow(data$board_dig_in_raw), byrow=TRUE)
    }
    
    if (header$num_board_dig_out_channels > 0) {
        data$board_dig_out_raw[1, indices$board_dig_out:(indices$board_dig_out + 128 - 1)] = matrix(readNUInt16(fid, 128), nrow=nrow(data$board_dig_out_raw), byrow=TRUE)
    }
               
    return(data)
}

# Define notch_filter function
notch_filter = function(input, f_sample, f_notch, bandwidth) {
    t_step = 1 / f_sample
    f_c = f_notch * t_step
    
    l = length(input)
    
    # Calculate IIR filter parameters
    d = exp(-2 * pi * (bandwidth / 2) * t_step)
    b = (1 + d * d) * cos(2 * pi * f_c)
    a0 = 1
    a1 = -b
    a2 = d * d
    a = (1 + d * d) / 2
    b0 = 1
    b1 = -2 * cos(2 * pi * f_c)
    b2 = 1
    
    output = c(input[1], input[2])
    
    # (If filtering a continuous data stream, change output[1] and output[2] to the previous final two values of out.)
    # Run filter
    for (k in 3:l) {
        output[k] = (a*b2*input[k-2] + a*b1*input[k-1] + a*b0*input[k] - a2*output[k-2] - a1*output[k-1])/a0
    }
    return(output)
}

# Define find_channel_in_group function
find_channel_in_group = function(channel_name, signal_group) {
    if (length(signal_group) > 0) {
        for (i in 1:length(signal_group)) {
            if (signal_group[[i]]$custom_channel_name == channel_name) {
                return(list("channel_found" = TRUE, "channel_index" = i))
            }
        }
    }
    return(list("channel_found" = FALSE, "channel_index" = 0))
}


# Define find_channel_in_header function
find_channel_in_header = function(channel_name, header) {
    # Look through all present signal groups
    
    # 1. Look through amplifier_channels
    returnList = find_channel_in_group(channel_name, header$amplifier_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "amplifier_channels", "signal_index" = channel_index))
    }
    
    # 2. Look through dc_amplifier_channels
    if (header$dc_amplifier_data_saved > 0) {
        returnList = find_channel_in_group(channel_name, header$dc_amplifier_channels)
        channel_found = returnList$channel_found
        channel_index = returnList$channel_index
        if (channel_found) {
            return(list("channel_found" = TRUE, "signal_type" = "dc_amplifier_channels", "signal_index" = channel_index))
        }
    }
    
    # 3. Look through stim_channels
    returnList = find_channel_in_group(channel_name, header$stim_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "stim_channels", "signal_index" = channel_index))
    }
    
    # 3.1 Look through amp_settle_channels
    returnList = find_channel_in_group(channel_name, header$amp_settle_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "amp_settle_channels", "signal_index" = channel_index))
    }
    
    # 3.2 Look through charge_recovery_channels
    returnList = find_channel_in_group(channel_name, header$charge_recovery_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "charge_recovery_channels", "signal_index" = channel_index))
    }
    
    # 3.3 Look through compliance_limit_channels
    returnList = find_channel_in_group(channel_name, header$compliance_limit_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "compliance_limit_channels", "signal_index" = channel_index))
    }
    
    # 4. Look through board_adc_channels
    returnList = find_channel_in_group(channel_name, header$board_adc_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "board_adc_channels", "signal_index" = channel_index))
    }
    
    # 5. Look through board_dac_channels
    returnList = find_channel_in_group(channel_name, header$board_dac_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "board_dac_channels", "signal_index" = channel_index))
    }
    
    # 6. Look through board_dig_in_channels
    returnList = find_channel_in_group(channel_name, header$board_dig_in_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "board_dig_in_channels", "signal_index" = channel_index))
    }
    
    # 7. Look through board_dig_out_channels
    returnList = find_channel_in_group(channel_name, header$board_dig_out_channels)
    channel_found = returnList$channel_found
    channel_index = returnList$channel_index
    if (channel_found) {
        return(list("channel_found" = TRUE, "signal_type" = "board_dig_out_channels", "signal_index" = channel_index))
    }
    
    return(list("channel_found" = FALSE, "signal_type" = "", "signal_index" = 0))
}

# Define plot_channel function
plot_channel = function(channel_name, result) {
    # Find channel that correpsonds to this name
    returnList = find_channel_in_header(channel_name, result)
    channel_found = returnList$channel_found
    signal_type = returnList$signal_type
    signal_index = returnList$signal_index
    
    # Plot this channel
    if (channel_found) {
        
        if (signal_type == "amplifier_channels") {
            y_label = "Voltage (microVolts)"
            data_vector = result$amplifier_data[signal_index, 1:ncol(result$amplifier_data)]
        } else if (signal_type == "dc_amplifier_channels") {
            y_label = "Voltage (Volts)"
            data_vector = result$dc_amplifier_data[signal_index, 1:ncol(result$dc_amplifier_data)]
        } else if (signal_type == "stim_channels") {
            y_label = "Current (microAmps)"
            data_vector = result$stim_data[signal_index, 1:ncol(result$stim_data)]
        } else if (signal_type == "amp_settle_channels") {
            y_label = "Amp Settle Events (High or Low)"
            data_vector = result$amp_settle_data[signal_index, 1:ncol(result$amp_settle_data)]
        } else if (signal_type == "charge_recovery_channels") {
            y_label = "Charge Recovery Events (High or Low)"
            data_vector = result$charge_recovery_data[signal_index, 1:ncol(result$charge_recovery_data)]
        } else if (signal_type == "compliance_limit_channels") {
            y_label = "Compliance Limit Events (High or Low)"
            data_vector = result$compliance_limit_data[signal_index, 1:ncol(result$compliance_limit_data)]
        } else if (signal_type == "board_adc_channels") {
            y_label = "Voltage (Volts)"
            data_vector = result$board_adc_data[signal_index, 1:ncol(result$board_adc_data)]
        } else if (signal_type == "board_dac_channels") {
            y_label = "Voltage (Volts)"
            data_vector = result$board_dac_data[signal_index, 1:ncol(result$board_dac_data)]
        } else if (signal_type == "board_dig_in_channels") {
            y_label = "Digital In Events (High or Low)"
            data_vector = result$board_dig_in_data[signal_index, 1:ncol(result$board_dig_in_data)]
        } else if (signal_type == "board_dig_out_channels") {
            y_label = "Digital Out Events (High or Low)"
            data_vector = result$board_dig_out_data[signal_index, 1:ncol(result$board_dig_out_data)]
        } else {
            stop(paste("Plotting not possible; signal type", signal_type, " not found", sep = ""))
        }
        
        plot(result$t, data_vector, type="l", main=channel_name, xlab="Time (s)", ylab=y_label)
    }
    
    else {
        stop(paste("Plotting not possible; channel ", channel_name, " not found", sep = ""))
    }
}


# Define readNInt32 function
readNInt32 = function(fid, nElements) {
    return(replicate(nElements, readInt32(fid)))
}

# Define readNUInt32 function
readNUInt32 = function(fid, nElements) {
    return(replicate(nElements, readUInt32(fid)))
}
                                    
# Define readNUInt16 function
readNUInt16 = function(fid, nElements) {
    return(replicate(nElements, readUInt16(fid)))
}

# Define readFloat32 function
readFloat32 = function(fid) {
    return(readBin(fid, numeric(), n=1, size="4", endian="little"))
}

# Define readUInt16 function
readUInt16 = function(fid) {
    return(sum(2^.subset(0:15, as.logical(rawToBits(readBin(fid, raw(), n=2, endian="little"))))))
}

# Define readUInt32 function
readUInt32 = function(fid) {
    return(sum(2^.subset(0:31, as.logical(rawToBits(readBin(fid, raw(), n=4, endian="little"))))))
}


# Define readInt16 function
readInt16 = function(fid) {
    return(readBin(fid, integer(), n=1, size="2", endian="little"))
}
            
# Define readInt32 function
readInt32 = function(fid) {
    return(readBin(fid, integer(), n=1, size="4", endian="little"))
}

# Define readQString function
readQString = function(fid) {
    # Read Qt style String. The first 32-bit unsigned number indicates the length of the string (in bytes).
    # If this number equals 0xffffffff, the string is null
    a = ""
    length = readUInt32(fid)
    if (length == 0xffffffff) {
        return("")
    }
    
    # Convert length from bytes to 16-bit Unicode words
    length = length / 2
    if (length > 0) {
        for (i in 1:length) {
            thisInt = readUInt16(fid)
            mode(thisInt) = "raw"
            thisChar = sapply(thisInt, rawToChar)
            a = paste(a, thisChar, sep = "")
        }
    }
    return(a)
}

# Define plural function
plural = function(n) {
    # Utility function to optionally pluralize words based on the value of n
    return(if (n == 1) "" else "s")
}
