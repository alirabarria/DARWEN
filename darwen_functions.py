import numpy as np
import pandas as pd
import subprocess
import os
import shutil
from pathlib import Path
import gc
from sklearn.metrics import r2_score, mean_absolute_percentage_error

file_path = os.path.dirname(os.path.realpath(__file__))

  

def list_archives(strAddress, format=".dat", key='', sort_data=False):
    """
    Lists all files in a given directory that match a specified format and contain a given key in their names.

    Parameters:
    - strAddress (str): The directory path where files are located.
    - format (str): The file extension to match. Defaults to ".dat".
    - key (str): An optional substring that should be present in the file names. Defaults to an empty string.

    Returns:
    - list: A sorted list of paths (as strings) to the files that match the given criteria.
    """
    directory = Path(strAddress).resolve()
    # Find matching files and construct full paths
    database = [str(file) for file in directory.iterdir() if (key in str(file) or key == '') and (
        format in str(file) or format=='')]


    # Sort the list of file paths
    if sort_data:
        return np.sort(database)
    else:
        return database
    


def correct_scientific_notation(input_str):
    from re import compile
    # Define a regex pattern to match the incorrect format
    pattern = compile(r"(\d+\.\d+)([+-])(\d+)")
    
    # Function to replace the incorrect format with the correct one
    def replace_func(match):
        base = match.group(1)
        sign = match.group(2)
        exponent = match.group(3)
        return f"{base}E{sign}{exponent}"
    
    # Apply the regex replacement
    corrected_str = pattern.sub(replace_func, input_str)
    
    return corrected_str


def out_to_dfs(dir_, error_purposes=False, file_key='_out', header_keyword='!', return_only_mols=True,
               big_err_def=404*1e5):
    """
    Reads simulation output files from a specified directory and extracts data into pandas DataFrames.
    
    Parameters:
    - dir_ (str): The directory from which to read the output files.
    - error_purposes (bool): Flag to indicate if the function is being used for error analysis.
    - file_key (str): Key to identify relevant files.
    - header_keyword (str): Character or string that identifies the header line in the files.
    - return_only_mols (bool): If True, returns only molecular data, excluding other simulation parameters.
    
    Returns:
    - tuple: Depending on 'error_purposes', returns a DataFrame (or list of DataFrames) and error indicators.
    """

    files = np.sort(os.listdir(dir_))
    matching_files = [file for file in files if file_key in file]
    

    dfs = []  # DataFrames extracted from files
    times, irate, rvirate, rtype, errs = ([] for i in range(5))
    cons_err = False

    if len(matching_files)==0:
        if error_purposes:
            return None, errs, 1e5
        else:
            return dfs, times, errs, irate, rvirate, rtype

    for file in matching_files:
        with open(os.path.join(dir_, file), 'r') as f:
            lines = f.readlines()
            # Find header line index based on the header keyword
            header_i = int(np.sum([header_keyword in line for i, line in enumerate(lines)])-1)

            # Extract error and conservation error
            err_line = lines[16].replace(' ', '')
            err_str = err_line.split('=')[1]
            try:
                err = float(err_str)
                
            except ValueError:
                err = float(correct_scientific_notation(err_str))
            except Exception as e:
                print (repr(e))
                err = big_err_def
            errs.append(err)

            cons_line = lines[20].replace(' ', '').replace('\n', '')
            if 'Conservation' in cons_line:
                cons_err = cons_line.split(':')[1] == '1' or cons_err

            if not error_purposes:
                # Process time and reaction rates if not for error purposes
                time_line = lines[15].replace(' ', '').replace('s', '')
                time = float(time_line.split('=')[1])
                times.append(time)

                if ('_i' in file and not '_i0' in file):
                    irate_line = lines[21].replace(' ', '')
                    irate.append(float(irate_line.split(':')[1]))
                    reac_type_line = lines[24].replace(' ', '')
                    rvirate_line = lines[25].replace(' ', '')
                    rtype.append(int(reac_type_line.split(':')[1]))
                    if ('reaction' in rvirate_line and rtype[-1] < 0):
                        rvirate.append(float(rvirate_line.split(':')[1]))
                    else:
                        rvirate.append(0)

                elif ('_i' in file and '_i0' in file):
                    irate.append(0)
                    rvirate.append(0)
                    rtype.append(0)

            df = pd.read_csv(os.path.join(dir_, file), header=header_i, sep='\s+')
            cols = np.append(df.columns[1:], 'nan')
            df.columns = cols
            if return_only_mols:
                # Drop columns not related to molecular data if specified
                df = df.drop(['height[km]', 'pressure[bar]', 'temperature[K]', 'nan'], axis=1)
            dfs.append(df)

    if error_purposes:
        return dfs[np.argmin(errs)], cons_err, np.min(errs)
    else:
        return dfs, np.array(times), np.array(errs), np.array(irate), np.array(rvirate), np.array(rtype)

    
    
    
def calculate_error(df_to_test, df_to_match, imp_molecs=[], mode='max_percent', big_loss_def=1e5, punish_by_mols=True, punish_cons_err=False,
                    slice_cols=False, major_threshold=1e-20, verbose=False, 
                    ignore_cols = ['height[km]', 'pressure[bar]','temperature[K]', 'MMW[g]', 'nan']):
    """
    Calculates the error between two dataframes based on specified metrics and important molecules.

    Parameters:
    - df_to_test: DataFrame containing test data.
    - df_to_match: DataFrame containing data to match against.
    - imp_molecs (list): List of important molecules to focus on. If empty, all molecules are considered.
    - mode (str): Error calculation mode. Supported modes: 'mre', 'mape', 'r2', 'max_percent'.
    - big_loss_def (float): Default high error value for penalization purposes.
    - punish_by_mols (bool): Flag to enable punishment based on missing important molecules.
    - punish_cons_err (bool): Flag to enable punishment if conservation errors are detected.
    - return_all_mols (bool): Flag to return a dictionary of errors of all molecules.
    - slice_cols (bool): Flag to slice columns before calculating error.

    Returns:
    - err (float): Calculated error based on the selected mode.
    - max_err (float): Maximum absolute error found.
    - mol_max_err (str): Molecule with the maximum error.
    - max_layer (int) : Layer in which the max error was found
    - mols_errors (dict): A dictionary of individual molecule errors
    """
    if df_to_test is None or punish_cons_err:
        return big_loss_def, big_loss_def, '', -1, {}
    
    mols = df_to_test.columns[3:-1] if slice_cols else df_to_test.columns

    overall_err, max_err, mol_max_err, max_layer = 0, 0, '', 0
    mol_errors = {}  # To store individual molecule errors

    for m in mols:
        if m in ignore_cols:
            continue
        is_in_test = m in df_to_test.columns
        is_in_true = m in df_to_match.columns
        if verbose:
            print (m, is_in_test, is_in_true)
        in_imp_mols = m in imp_molecs or len(imp_molecs)==0
        if is_in_test and is_in_true:  
            y_true = np.power(10, df_to_match[m].values)  # Assuming log-transformed values
            y_test = np.power(10, df_to_test[m].values)
            diff = np.abs(y_true - y_test)
            max_diff = np.max(diff)
            rel_loss = diff/y_true
            major_array = rel_loss[y_true > major_threshold]
            if len(major_array)>0:
                major_loss = np.max(major_array)
            else:
                major_loss = 0

            if max_diff > max_err:
                max_err, mol_max_err, max_layer = max_diff, m, np.argmax(diff)

            # Calculate and store individual molecule errors
            if mode == 'mre':
                mol_error = np.mean(rel_loss)
                mol_errors[m] = [mol_error, max_diff, major_loss]
                if in_imp_mols:
                    overall_err+=mol_error

            elif mode == 'max_percent':
                mol_error = np.max(rel_loss)
                mol_errors[m] = [mol_error, max_diff, major_loss]
                if in_imp_mols:
                    overall_err = np.max([overall_err, mol_error])
           
      
        elif punish_by_mols and in_imp_mols:
            print ('Problem in mol:', m)
            overall_err = 237*big_loss_def
            continue
        else:
            overall_err = 235*big_loss_def
            max_err = big_loss_def
            mol_max_err = ''

    
    if mode == 'mre':
        if len(imp_molecs)==0:
            overall_err = overall_err / len(mols)
        else:
            overall_err = overall_err / len(imp_molecs)

    return overall_err, max_err, mol_max_err, max_layer, mol_errors

def from_df_to_array(elem):
    if isinstance(elem, str):
        return eval(elem)
    else:
        return elem
    
def list_directories_with_string(path, search_string):
    matching_directories = []

    # Check if the path exists
    if not os.path.exists(path):
        print(f"Path '{path}' does not exist.")
        return matching_directories

    # Get a list of all items in the directory
    items = os.listdir(path)

    # Iterate over each item and check if it is a directory and contains the search string
    for item in items:
        item_path = os.path.join(path, item)
        if os.path.isdir(item_path) and search_string in item:
            matching_directories.append(item_path)

    return matching_directories





def move_files_with_string(folder_path, target_folder, string_to_match, delete_previous=False,
                          verbose=False, copy=False):
    """
    Moves files that contain a specific string in their names from one folder to another.

    Parameters:
    - folder_path (str): The source folder path where files are located.
    - target_folder (str): The target folder path where files will be moved.
    - string_to_match (str): The substring to look for in file names.
    - delete_previous (bool): If True, deletes all files in the target folder before moving new files.
    - verbose (bool): If True, prints out the actions being taken.

    Returns:
    None: Files are moved from the source to the target folder, and no value is returned.
    """
    # Check if the target folder exists, create it if not
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
        if verbose:print(f'Created target folder: {target_folder}')
    elif delete_previous:
        # If delete_previous is True, delete the contents of the target folder
        if verbose:print(f'Deleting previous contents of {target_folder}')
        shutil.rmtree(target_folder)
        os.makedirs(target_folder)
    
    # Use the find command in bash to locate files containing the specified string
    bash_command = f'find {folder_path} -type f -name "*{string_to_match}*"'
    env = os.environ.copy()
    # Run the bash command and capture the output
    result = subprocess.run(bash_command, shell=True, capture_output=True, text=True, env=env)

    # Extract the file paths from the output
    file_paths = result.stdout.split('\n')

    # Remove any empty strings from the list
    file_paths = [path for path in file_paths if path]

    if file_paths:
        # Move the files to the target folder
        for file_path in file_paths:
            move_command = f'mv {file_path} {target_folder}' if not copy else f'cp {file_path} {target_folder}'
            subprocess.run(move_command, shell=True)
        if verbose:print(f'Moved {len(file_paths)} files to {target_folder}')
    elif verbose:
        print('No files found with the specified string.')

def extract_molecs(line, mode='array', link='->', return_arrehnius=False, int_max=28):
    """
    Extracts molecules from a reaction line, supporting different output formats.

    Parameters:
    - line (str): The reaction line from which to extract molecules.
    - mode (str): Determines the output format ('array' for lists of reactants and products, 'text' for a formatted string).
    - link (str): The linking symbol between reactants and products in 'text' mode. Default is '->'.
    - check_line (int): Additional characters to check beyond the reaction, used in 'array' mode with a nonzero value.

    Returns:
    Depending on 'mode', either lists of reactants and products, possibly including additional checked characters, or a formatted string representation of the reaction.
    """
    try:
        sep = line.split('=')
        reacs = sep[0].strip().split('+')
        reacs = [r.strip() for r in reacs]
        prods = sep[1][:int_max].strip().split('+')
        prods = [p.strip() for p in prods]
    except Exception as e:
        print (line)
    
    if mode == 'array' and not return_arrehnius:
        return reacs, prods
    
    elif mode=='array' and return_arrehnius:
        arr = line.split(':')
        splitted = arr[0].split()
        rtype = int(splitted[-1])
        gamma = np.float64(splitted[-2])
        beta = np.float64(splitted[-3])
        alpha = np.float64(splitted[-4])
        return reacs, prods, rtype, alpha, beta, gamma

    
    elif mode == 'text':
        to_return = ""
        first = True
        for rct in reacs:
            if first:
                to_return+= rct
                first=False
            else:
                to_return+=' + '+rct
        to_return+=' '+link+' '
        first = True
        for pr in prods:
            if first:
                to_return+= pr
                first=False
            else:
                to_return+=' + '+pr
        return to_return

def network_to_df(reacs_full_fname='1D_PAC/venot2020.reac', force_id=False,
           check_tol = 1e-4, signal_wrong_ids=True, id_col='id', return_df=True,
           search_id = False, other_reacs_df=[], original_data_df=[], find_Fg=True):
    reacs_dict = {}
    reacs_file = open(reacs_full_fname, 'r')
    reacs_lines = reacs_file.readlines()
    to_dict = ''
    c = 0
    df = pd.DataFrame()
    for i in range(len(reacs_lines)): 
        not_spaces = reacs_lines[i].replace(' ', '')
        if len(not_spaces.replace('\n', '')) == 0 or not_spaces[0] == '!':
            continue
        to_dict += reacs_lines[i]
        if i < len(reacs_lines) - 1:
            next_line = reacs_lines[i+1]
        else:
            next_line = ''

        if len(next_line[:4].replace(' ', '')) != 0 or i + 1 == len(reacs_lines):
            reacs_dict[c] = to_dict
            to_dict = ''
            c+=1

        if return_df and to_dict == '':
            ind = c - 1
            reacs, prods, rtype, alpha, beta, gamma = extract_molecs(reacs_dict[c - 1], return_arrehnius=True)
            df.at[ind, 'reacs'], df.at[ind, 'prods'], df.at[ind, 'rtype'] = str(reacs), str(prods), rtype
            df.at[ind, 'molecs'] = np.append(reacs, prods)
            df.at[ind, 'alpha'], df.at[ind, 'beta'], df.at[ind, 'gamma'] = alpha, beta, gamma

            if abs(rtype) in [800, 810, 811, 812, 900, 910, 911]:
                df.at[ind, 'k_inf'] = str(reacs_dict[c - 1].split('\n')[1].split()[:3])
                df.at[ind, 'Troe'] = str(reacs_dict[c - 1].split('\n')[2].split()[:4])
                feff = {}
                for line in reacs_dict[c - 1].split('\n')[3:]:
                    if line == '':
                        continue
                    feff[line.split()[0]] = line.split()[1]
                df.at[ind, 'feff'] = str(feff)
            
            elif abs(rtype) in [801, 901]:
                feff = {}
                for line in reacs_dict[c - 1].split('\n')[2:]:
                    if line == '':
                        continue
                    feff[line.split()[0]] = line.split()[1]
                df.at[ind, 'feff'] = str(feff)
            
            elif (30 < abs(rtype) and abs(rtype) < 45) or abs(rtype)==908:
                prefs = [line.split()[-1] if len(line.replace(' ', '')) != 0 else '' for line in reacs_dict[c - 1].split('\n')]
                if prefs[-1] == '':
                    prefs = prefs[:-1]
                alphas, betas, gammas = [alpha], [beta], [gamma]
                for line in reacs_dict[c - 1].split('\n')[1:]:
                    if line == '':
                        continue
                    alphas.append(line.split()[0])
                    betas.append(line.split()[1])
                    gammas.append(line.split()[2])
                df.at[ind, 'prefs'] = str(prefs)
                df.at[ind, 'alphas'] = str(alphas)
                df.at[ind, 'betas'] = str(betas)
                df.at[ind, 'gammas'] = str(gammas)

            
            found_id = False
            if search_id and len(other_reacs_df) > 0:
                df_candidates = other_reacs_df[(other_reacs_df.reacs==str(reacs)) & (other_reacs_df.prods==str(prods))]
                if len(df_candidates) > 0:
                    ind = df_candidates.index[0]
                    alpha_row = df_candidates.iloc[0]['alpha']
                    if abs(alpha - float(alpha_row)) < check_tol or force_id:
                        other_reacs_df.at[ind, id_col] = c
                        found_id = True
                    elif signal_wrong_ids:
                        print ('signaling:', ind, c - 1, alpha_row, alpha)

            if len(original_data_df) > 0 and find_Fg:
                try:
                    df.at[ind, 'F'], df.at[ind, 'g'] = original_data_df.at[ind, 'F'], original_data_df.at[ind, 'g']
                         
                except Exception as e:
                    print ("Couldn't find F and g for", ind)
            
    if return_df and search_id:
        return reacs_dict, df, other_reacs_df
    elif return_df:
        return reacs_dict, df
    else:
        return reacs_dict


"""
def add_id(reacs_in_venot=[], reacs_full_fname='1D_PAC/venot2020.reac', force_id=False,
           check_tol = 1e-4, signal_wrong_ids=True, id_col='id', mid_id = 837):
    reacs_from_file = {}
    reacs_file = open(reacs_full_fname, 'r')
    reacs_lines = reacs_file.readlines()
    i = 0
    c = 0
    while i < len(reacs_lines):
        if i < mid_id:
            reacs_from_file[c] = reacs_lines[i]
            c+=1
            i+=1
            
        else:
            to_dict = reacs_lines[i]
            i+=1
            if i < len(reacs_lines):
                aux = reacs_lines[i]
                while len(aux[:4].replace(' ', '')) == 0:
                    to_dict+=aux
                    i+=1
                    aux = reacs_lines[i]

                reacs_from_file[c] = to_dict       
                c+=1
        if len(reacs_in_venot) > 0:
            reacs, prods, rtype, alpha, beta, gamma = extract_molecs(reacs_from_file[c-1], return_arrehnius=True)
            venot_candidates = reacs_in_venot[(reacs_in_venot.reacs==str(reacs)) & (reacs_in_venot.prods==str(prods))]
            if len(venot_candidates) > 0:
                ind = venot_candidates.index[0]
                alpha_row = venot_candidates.iloc[0]['alpha']
                if abs(alpha - float(alpha_row)) < check_tol or force_id:
                    reacs_in_venot.at[ind, id_col] = c - 1
                elif signal_wrong_ids:
                    print ('signaling:', ind, c - 1, alpha_row, alpha)
    return reacs_in_venot, reacs_from_file
"""

def add_id(reacs_in_venot=[], reacs_full_fname='1D_PAC/venot2020.reac', force_id=False,
           check_tol = 1e-4, signal_wrong_ids=True, id_col='id', mid_id = 837):
    reacs_from_file = {}
    reacs_file = open(reacs_full_fname, 'r')
    reacs_lines = reacs_file.readlines()
    to_dict = ''
    c = 0
    for i in range(len(reacs_lines)-1): 
        to_dict += reacs_lines[i]
        next_line = reacs_lines[i+1]
        if len(next_line[:4].replace(' ', '')) != 0 or i + 1 == len(reacs_lines):
            reacs_from_file[c] = to_dict
            c+=1
            to_dict = ''

        if len(reacs_in_venot) > 0:
            reacs, prods, rtype, alpha, beta, gamma = extract_molecs(reacs_from_file[c-1], return_arrehnius=True)
            venot_candidates = reacs_in_venot[(reacs_in_venot.reacs==str(reacs)) & (reacs_in_venot.prods==str(prods))]
            if len(venot_candidates) > 0:
                ind = venot_candidates.index[0]
                alpha_row = venot_candidates.iloc[0]['alpha']
                if abs(alpha - float(alpha_row)) < check_tol or force_id:
                    reacs_in_venot.at[ind, id_col] = c - 1
                elif signal_wrong_ids:
                    print ('signaling:', ind, c - 1, alpha_row, alpha)
    return reacs_in_venot, reacs_from_file

"""
def run_pac(pop_str, silent=False, hide_output=False, out_file=None, out_path=Path(), max_processes=6, command='./'):

    if not silent:
            print('.........................................................................................................', file=out_file)
            print('                                       Running', pop_str, file=out_file)
            print('.........................................................................................................', file=out_file)
            if hide_output: out_file.flush()

    # Running PAC
    output_file = (out_path / Path(f'output_{pop_str}.txt')).resolve()
    # Ensure the parent directory exists
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # Create the file
    output_file.touch(exist_ok=True)

    command = f'{command}multipac.sh {pop_str} {max_processes} > {output_file} 2>&1'
    gc.collect()
    with open(output_file, 'w') as output_file:
        process = subprocess.Popen(command, shell=True, stdout=output_file, stderr=subprocess.STDOUT, 
                                       cwd=out_path)
        process.wait()  # Wait for the process to finish
"""


def run_pac(pop_str, silent=False, hide_output=False, out_file=None, out_path=Path(), max_processes=6, command='./'):
    if not silent:
        print('.........................................................................................................', file=out_file)
        print('                                       Running', pop_str, file=out_file)
        print('.........................................................................................................', file=out_file)
        if hide_output:
            out_file.flush()

    # Prepare output file path
    output_file_path = (out_path / f'output_{pop_str}.txt').resolve()
    output_file_path.parent.mkdir(parents=True, exist_ok=True)

    # Build command as list (safe, no shell=True)
    command_list = [command + 'multipac.sh', pop_str, str(max_processes)]

    gc.collect()

    with open(output_file_path, 'w') as output_file:
        process = subprocess.Popen(
            command_list,
            stdout=output_file,
            stderr=subprocess.STDOUT,
            cwd=out_path
        )
        process.wait()


def format_line(*args, total_length=32):
    """
    Formats a line to include given numerical values without changing the total length of the line.

    Parameters:
    args (tuple): Numerical values to be formatted into the string.
    total_length (int): The fixed length of the string.

    Returns:
    str: A formatted string of fixed length.
    """
    # Convert all arguments to strings
    str_args = [str(arg) for arg in args]
    
    # Initial check if the total length of numbers and minimum spaces exceeds the limit
    if sum(len(num) for num in str_args) + len(str_args) - 1 > total_length:
        raise ValueError("The numbers and required spaces exceed the total length.")

    # Calculate the total number of spaces available after placing the numbers
    remaining_spaces = total_length - sum(len(num) for num in str_args)

    # Determine the minimum number of spaces between elements
    if len(str_args) > 1:
        min_spaces_between = remaining_spaces // (len(str_args) - 1)
    else:
        min_spaces_between = remaining_spaces

    # Construct the line with evenly distributed spaces
    formatted_line = ''
    for i, num in enumerate(str_args):
        if i < len(str_args) - 1:
            # Add the number and the spaces after it
            formatted_line += num + ' ' * min_spaces_between
            remaining_spaces -= min_spaces_between
        else:
            # Add the last number and any remaining spaces
            formatted_line += num

    # If there are any leftover spaces, add them to the end
    formatted_line += ' ' * remaining_spaces

    return formatted_line

from datetime import datetime, timedelta

def get_yesterday(input_date):
    # Convert the string date to a datetime object
    date = datetime.strptime(input_date, "%Y-%m-%d")
    # Subtract one day
    yesterday = date - timedelta(days=1)
    # Return the date in the same string format
    return yesterday.strftime("%Y-%m-%d")

def mol_to_latex(molec):
    if molec == 'H2O':
        return '\mathrm{H}_2 \mathrm{O}'
    elif molec == 'CH4':
        return '\mathrm{CH}_4'
    elif molec == 'CO2':
        return '\mathrm{CO}_2'
    elif molec=='C2H2':
        return '\mathrm{C}_2 \mathrm{H}_2'
    elif molec == 'NH3':
        return '\mathrm{NH}_3'
    else:
        return molec
    
from decimal import Decimal, getcontext
def arrhenius(t, alpha, beta, gamma, prec=50):
    getcontext().prec=prec
    t, alpha, beta, gamma =Decimal(t), Decimal(alpha), Decimal(
        beta), Decimal(gamma)
    return Decimal(alpha*(t/300)**beta*np.exp(-gamma/t))

def get_rate_constant(t, p, rtype, alpha=1, beta=0, gamma=0, prefs=[], feff={},
                      extra_params=[], inf_params=[], abuns={}, kboltz=1.3806504e-16,
                      eps=1e-300, prec=30):
    getcontext().prec = prec
    abs_rtype = abs(rtype)
    p = p *1e6
    density = Decimal(p/(kboltz*t))
    t = Decimal(t)
    p = Decimal(p)

    if abs_rtype>30 and abs_rtype < 40:
        nrefpr = int(abs_rtype - 30)
        if len(alpha) == 1:
            print ('Plog reaction not well handled')
        if p < prefs[0]: 
            k_rate = arrhenius(t, alpha[0], beta[0], gamma[0], prec=prec)
        elif p >=prefs[-1]:
            k_rate = arrhenius(t, alpha[-1], beta[-1], gamma[-1], prec=prec)
        else:
            for ref in range(1, nrefpr-1):
                p1 = prefs[ref]
                p2 = prefs[ref+1]
                if p <= p2:
                    klow =  arrhenius(t, alpha[ref], beta[ref], gamma[ref], prec=prec)
                    khigh = arrhenius(t, alpha[ref+1], beta[ref+1], gamma[ref+1], prec=prec)

                    logk = np.log(klow) + (np.log(khigh) - np.log(klow)
                                           )*(np.log(p) - np.log(p1)
                                              )/(np.log(p2) - np.log(p1))
                    k_rate = np.exp(logk)
    elif abs_rtype in [800, 810, 811, 812, 900, 910, 911]:
        k0 = arrhenius(t, alpha, beta, gamma, prec=prec)
        kinf = arrhenius(t, inf_params[0], inf_params[1], inf_params[2], prec=prec)
        m=0
        if abs_rtype == 811:
            m=abuns['H2O']*density
        elif abs_rtype == 812:
            m=abuns['CO2']*density
        elif abs_rtype == 911:
            m=abuns['CH4']*density
        else:
            for mol in abuns.keys():
                if mol in feff.keys():
                    f_eff = Decimal(feff[mol])
                    abun = Decimal(abuns[mol])
                    m=m+f_eff*abun*density

        if abs_rtype in [800, 811, 812, 900, 911]: #Troe formulation
            a, t3, t1, t2, N = extra_params
            a, t3, t1, t2 = Decimal(a), Decimal(t3 + eps), Decimal(t1+eps), Decimal(t2)
            fcent = (1 - a)*np.exp(-t/t3) + a*np.exp(-t/t1) + np.exp(-t2/t)
            fcent = np.max([Decimal(eps), fcent])
            c=-Decimal(0.4)-Decimal(0.67)*np.log10(fcent)
            n=Decimal(0.75)-Decimal(1.27)*np.log10(fcent)
            d=Decimal(0.14)
            if m > 0:
                fk = np.log10(fcent)/(
                    1+((np.log10(k0*m/kinf)+c)/n-d*np.log10(k0*m/kinf)+c)**Decimal(2.))
                fk = np.power(10, fk)
            else:
                fk=0
        elif abs_rtype in [810, 910]:
            a, b, c, d, e, N = extra_params
            x=1/(1+(np.log10(k0*m/kinf))**2.)
            fk=d*t**e*(a*np.exp(-b/t)+np.exp(-t/c))**x
        k_rate = k0*m/(1+k0*m/kinf)*fk

    elif abs_rtype in [801, 901]:
        k0=arrhenius(t, alpha, beta, gamma)
        m = 0
        for mol in abuns.keys():
            if mol in feff.keys():
                m=m+feff[mol]*abuns[mol]*density
        k_rate = k0*m
    else:
        k_rate = arrhenius(t, alpha, beta, gamma)

    return k_rate

def get_k(reacs_df, ind, t, p, abuns={}):
    data = reacs_df.at[ind, 'data']
    rtype = reacs_df.at[ind, 'rtype']
    abs_rtype = abs(rtype)
    unc_type = ''
    feff = {}
    alpha = reacs_df.at[ind, 'alpha']
    beta = reacs_df.at[ind, 'beta']
    gamma = reacs_df.at[ind, 'gamma']
    
    if (30 < abs_rtype and abs_rtype < 45) or abs(rtype)==908: #Plog
        prefs = []
        array_data = data.split('\n')[:-1]
        unc_type = array_data[0].split(':')[1].strip()
        prefs.append(float(array_data[0].split(':')[2].strip()))
        alpha = [alpha]
        beta = [beta]
        gamma = [gamma]
        for ref in range(1, len(array_data)):
            aux_data = array_data[ref]
            pointy = aux_data.split(':')
            prefs.append(float(pointy[1].strip()))
            spacey = pointy[0].split()
            alpha.append(float(spacey[0]))
            beta.append(float(spacey[1]))
            gamma.append(float(spacey[2]))
        k_rate = get_rate_constant(t, p, rtype=rtype, alpha=alpha, beta=beta, gamma=gamma, prefs=prefs)

    elif abs_rtype in [800, 810, 811, 812, 900, 910, 911]: #Troe params
        feff = {}
        array_data = data.split('\n')[:-1]
        unc_type = array_data[0].split(':')[1].strip()
        aux_params = array_data[1].split(':')
        inf_params = aux_params[0].split()

        aux_params = array_data[2].split(':')
        troe_params = aux_params[0].split()

        inf_params = [float(aux) for aux in inf_params]
        troe_params = [float(aux) for aux in troe_params]
        for ref in range(3, len(array_data)):
            aux_params = array_data[ref].split(':')[0].split()
            feff[aux_params[0]] = float(aux_params[1])
        k_rate = get_rate_constant(t, p, rtype=rtype, abuns=abuns,
                                   alpha=alpha, beta=beta, gamma=gamma,
                                   feff=feff, inf_params=inf_params, 
                                   extra_params=troe_params)
    elif abs_rtype in [810, 910]:
        array_data = data.split('\n')[:-1]
        unc_type = array_data[0].split(':')[1].strip()
        for ref in range(2, len(array_data)):
            aux_params = array_data[ref].split(':')[0].split()
            feff[aux_params[0]] = float(aux_params[1])
        k_rate = get_rate_constant(t, p, rtype=rtype, feff=feff,
                                   abuns=abuns, alpha=alpha, beta=beta,
                                   gamma=gamma)
    else:
        unc_type = data.split(':')[1].strip()
        k_rate = arrhenius(t, alpha=alpha, beta=beta, gamma=gamma)


    return float(k_rate), unc_type

def calculate_rate(reacs_df, abun_df, ind, ind_height):
    row = reacs_df.iloc[ind]
    reacs = row['reacs']
    krate = row['krate'+str(ind_height)]
    nu = Decimal(krate)
    for reac in reacs:
        abun = abun_df.at[ind_height, reac]
        nu *= Decimal(10**abun)
    return float(nu)



def obtain_unc_factors(prefix = '', feff_keys = ['O2', 'CO', 'CO2', 'H2O', 'CH4', 'H2', 'C2H6', 
                                                     'Ar', 'N2', 'He', 'C2H4', 'benzene', 'toluene', 
                                                     'H2O2', 'N2O', 'O3P', 'NH3', 'N2H4', 'N2O4', 
                                                     'NO2', 'NO', 'N4S'],
                      columns=['reacs', 'prods', 'alpha', 'beta', 'gamma', 'Troe', 'SRI', 'k_inf', 'feff', 'prefs', 'alphas',
                                            'betas', 'gammas', 'F', 'g', 'F_inf', 'g_inf'] ):
    

    data_from_venot = pd.DataFrame(columns=columns)
    c = 0
    for i in range(1,14):
        direc = prefix+str(i)+'.dat'
        file_lines = open(direc, 'r').readlines()
        n = 0
        for line in file_lines:
            reacs = line[:56]
            prods = line[56:106]
            params = line[106:]
            data_from_venot.at[c, 'reacs']= str(reacs.split())
            data_from_venot.at[c, 'prods']=str(prods.split())
            params = params.split()

            if i in [1,2,8]: #Reversible, Irreversible bimolecular reacs, Reversible decomposition reactions
                alpha, beta, gamma, F, g = params
                
            elif i in [3,9]: #Reversible bimolecular with pressure dependence, Reversible decomposition reactions - coefficients with a pressure dependence

                if n < 3:
                    N = 5
                else:
                    N = 8
                if i==9:
                    N = 6
                prefs = params[:N]
                alphas = params[N:2*N]
                betas = params[2*N:3*N]
                gammas = params[3*N:4*N]
                F = params[4*N]
                g = params[4*N+1]
                alpha = alphas[0]
                beta=betas[0]
                gamma=gammas[0]
                data_from_venot.at[c, 'prefs'] = str(prefs)
                data_from_venot.at[c, 'alphas'] = str(alphas)
                data_from_venot.at[c, 'betas'] = str(betas)
                data_from_venot.at[c, 'gammas'] = str(gammas)
                n += 1

            elif i in [4,10,13]: #Reversible combination reactions, low pressure limit, Reversible decomposition reactions - low pressure limit, Reversible third-body assisted desexcitation reactions
                # Alpha_0   Beta_0   Gamma_0   F_0   g_0   Efficiencies
                feff = {}
                alpha, beta, gamma, F, g = params[:5]
                for j in range(len(params[5:])):
                    feff[feff_keys[j]] = params[5+j]
                data_from_venot.at[c, 'feff'] = str(feff)

            elif i in [5, 7, 11]: #Reversible combination reactions - Troe formalism, Reversible combination reactions with specific third-body - Troe formalism, Reversible decomposition reactions - Troe formalism
                # Alpha_0   Beta_0   Gamma_0   F_0   g_0   Alpha_inf   Beta_inf   Gamma_inf   F_inf  g_inf   A_Troe   B_Troe   C_Troe   D_Troe   Efficiencies
                alpha, beta, gamma, F, g = params[:5]
                data_from_venot.at[c, 'k_inf'] = params[5:8]
                data_from_venot.at[c, 'F_inf'] = params[8]
                data_from_venot.at[c, 'g_inf'] = params[9]
                data_from_venot.at[c, 'Troe'] = params[10:14]
                N = 14
                for j in range(len(params[N:])):
                    feff[feff_keys[j]] = params[N+j]
                data_from_venot.at[c, 'feff'] = str(feff)

            elif i in [6, 12]: #Reversible combination reactions - SRI formalism, Reversible decomposition reactions - SRI formalism
                # Alpha_0   Beta_0   Gamma_0   F_0   g_0   Alpha_inf   Beta_inf   Gamma_inf   F_inf  g_inf  A_SRI   B_SRI   C_SRI   D_SRI   E_SRI   Efficiencies
                alpha, beta, gamma, F, g = params[:5]
                data_from_venot.at[c, 'k_inf'] = params[5:8]
                data_from_venot.at[c, 'F_inf'] = params[8]
                data_from_venot.at[c, 'g_inf'] = params[9]
                data_from_venot.at[c, 'SRI'] = params[10:15]
                N = 15
                for j in range(len(params[N:])):
                    feff[feff_keys[j]] = params[N+j]
                data_from_venot.at[c, 'feff'] = str(feff)


            data_from_venot.at[c, 'alpha'] = alpha
            data_from_venot.at[c, 'beta'] = beta
            data_from_venot.at[c, 'gamma'] = gamma
            data_from_venot.at[c, 'F'] = F
            data_from_venot.at[c, 'g'] = g
            data_from_venot.at[c, 'num'] = i
            c+=1
    return data_from_venot








def add_Fg(reacs_df, data_from_venot, use_whole_net=True):
    for ind, row in reacs_df.iterrows():
        if use_whole_net:
            id_ = ind
        else:
            id_ = row['id']

        if id_ == 958:
            continue
        else:
            row_venot = data_from_venot.iloc[int(id_)]


        if str(row['reacs']) == str(row_venot['reacs']):
            F, g, F_inf, g_inf = row_venot['F'], row_venot['g'], row_venot['F_inf'], row_venot['g_inf']
            
        else:
            for ind_venot, row_venot in data_from_venot.iterrows():
                if row['reacs'] == eval(row_venot['reacs']
                        ) and row['prods']==eval(row_venot['prods']):
                    F, g, F_inf, g_inf = row_venot['F'
                                ], row_venot['g'], row_venot['F_inf'], row_venot['g_inf']
                else:
                    print ('problem with', ind, row['reacs'], row['prods'] )

        reacs_df.at[ind, 'F'] = F
        reacs_df.at[ind, 'g'] = g
        reacs_df.at[ind, 'F_inf'] = F_inf
        reacs_df.at[ind, 'g_inf'] = g_inf
    return reacs_df





def unc_factor(t, F0, g0, T0=300):
    return F0*np.exp(g0*np.abs(1/t - 1/T0))

def add_params(reacs_df, abun_df, add_krate=False):
    for ind in reacs_df.index:
        nu_max = 0
        unc_max = 0
        unc_min = 10
        for ind_abun, row_abun in abun_df.iterrows():
            t = row_abun['temperature[K]']
            p = row_abun['pressure[bar]']
            abuns = row_abun
            if add_krate:
                krate, unc_type = get_k(reacs_df, ind, t, p, 
                                    abuns=abuns)
                reacs_df.at[ind, 'krate'+str(ind_abun)] = krate
                nu = calculate_rate(reacs_df, abun_df, ind, ind_abun)
                reacs_df.at[ind, 'nu'+str(ind_abun)] = nu
                if nu > nu_max:
                    nu_max = nu
                

            F = float(reacs_df.at[ind, 'F'])
            g = float(reacs_df.at[ind, 'g'])
            sigma = unc_factor(t, F0=F, g0=g)

            reacs_df.at[ind, 'unc'+str(ind_abun)] = sigma
            

            if sigma > unc_max:
                unc_max = sigma
            
            if sigma < unc_min:
                unc_min = sigma
        if add_krate:
            reacs_df.at[ind,  'unc_type'] = unc_type
            reacs_df.at[ind, 'nu_max'] = nu_max
        reacs_df.at[ind, 'unc_max'] = unc_max
        reacs_df.at[ind, 'unc_min'] = unc_min
    return reacs_df




def initialize(network='veillet', fetuccini_mode=False, force_recalc=False,
                lira_required=False):
    
    global imp_molecs, abund_df, all_molecs, prefix, reacs_from_lira, reacs_in_lira, specs_from_file, specs_lines
    imp_molecs = ['H2O', 'CH4', 'CO', 'CO2', 'NH3', 'HCN', 'C2H2', 'OH']
    if lira_required:
        reacs_from_lira, reacs_in_lira = network_to_df(reacs_full_fname=os.path.join(file_path, '1D_PAC/lira_val2024.reac'), 
                                                        return_df=True)
    
    if network == 'veillet':
        
        global data_from_veillet, reacs_from_veillet, reacs_in_veillet
        print ('initializing veillet vars')
        all_molecs = np.load(os.path.join(file_path, '1D_PAC/mols_veillet.npy'))
        specs_file = open(os.path.join(file_path, '1D_PAC/veillet2023.spec'), 'r')
        data_from_veillet = pd.read_csv(os.path.join(file_path, '1D_PAC/data_from_veillet.csv'))
        reacs_from_veillet, reacs_in_veillet = network_to_df(reacs_full_fname=os.path.join(file_path, '1D_PAC/veillet2023.reac'),
                            force_id=True, id_col='veillet_id', original_data_df=data_from_veillet,
                            search_id = False, return_df=True)
                            
        prefix = os.path.join(file_path, '1D_PAC/VENOT2020/reactions_')

        if fetuccini_mode:
            vdf, cons_err, _ = out_to_dfs('1D_PAC/hd2094_veillet', error_purposes=True, file_key='_out', 
                                   return_only_mols=False)
            reacs_in_veillet = add_params(reacs_in_veillet, vdf, add_krate=False)
        
    else:
        global data_from_venot, reacs_from_venot, reacs_from_file, reacs_in_venot, reacs_in_full
        print ('initializing venot vars')
        prefix = '/Users/arturolira/phd_no_icloud/DARWEN/VENOT2020/reactions_'
        

        if os.path.exists(os.path.join(file_path, 'data_from_venot.csv')) and not force_recalc:
            data_from_venot = pd.read_csv(os.path.join(file_path, 'data_from_venot.csv'))
        else:
            data_from_venot = obtain_unc_factors(prefix = prefix)
            data_from_venot.to_csv(os.path.join(file_path, 'data_from_venot.csv'), index=False)
        all_molecs = np.genfromtxt(os.path.join(file_path,'1D_PAC/mols.txt'), dtype=None, encoding=None)
        specs_file = open(os.path.join(file_path, '1D_PAC/venot2020.spec'), 'r')
        reacs_from_venot, reacs_in_venot = network_to_df(reacs_full_fname=os.path.join(file_path, '1D_PAC/venot2020reduced.reac'), 
                                            return_df=True, find_Fg=False)
        reacs_from_file, reacs_in_full = network_to_df(reacs_full_fname=os.path.join(file_path, '1D_PAC/venot2020.reac'), return_df=True,
                    search_id = False, force_id=True, id_col='id', original_data_df=data_from_venot,
                    other_reacs_df=[], find_Fg=True)
        if fetuccini_mode:
            abun_df, cons_err, _ = out_to_dfs(os.path.join(file_path,'1D_PAC/hd2094_no_photo'), 
                    error_purposes=True, return_only_mols=False)
            if 'MMW[g]' in abun_df.columns:
                abun_df = abun_df.drop(['nan', 'MMW[g]'], axis=1)
            else:
                abun_df = abun_df.drop(['nan'], axis=1)
            reacs_in_full = add_params(reacs_in_full, abun_df, add_krate=False)

    specs_lines = specs_file.readlines()
    specs_from_file = {}
    for i in range(3,len(specs_lines)):
        sp = specs_lines[i].split(' ')[1]
        specs_from_file[sp] = specs_lines[i]

    specs_file.close()

   


def write_new_scheme(kid_reduc, reduc_molecs, fname='lira2023', verbose=False, loss=-1,
                    write_ids=False, imp_mols_flag=True, add_He_reac=False,
                    reacs_dict={}, specs_lines=[], specs_from_file={}):
    """
    Writes a new scheme into files based on reduced sets of reactions and molecules.

    Parameters:
    - kid_reduc (np.array): Array of reaction IDs to include in the new scheme.
    - reduc_molecs (np.array): Array of molecule IDs to include in the new scheme.
    - fname (str): Base filename for the output files (without extension).
    - verbose (bool): If True, prints additional information during the process.
    - loss (float): The loss value associated with this scheme, for documentation in the output.
    - write_ids (bool): If True, writes the reaction IDs into a separate file.
    - imp_mols_flag (bool): Indicates if the loss pertains only to important molecules.
    - add_He_reac (bool): If True, artificially adds a Helium reaction to the scheme.

    Returns:
    None: This function writes to files and does not return any value.
    """
    kid_reduc=np.unique(kid_reduc)
    reduc_molecs=np.unique(reduc_molecs)
    if write_ids:
        with open(fname + '.int', 'w') as int_file:
            np.savetxt(int_file, np.sort(kid_reduc))

    if verbose and loss ==-1:
        print ('Writing', len(kid_reduc), 'reacs and', len(reduc_molecs), 'species.')
    elif verbose:
        print ('Writing', len(kid_reduc), 'reacs and', len(reduc_molecs), 'species and MRE=', loss)
    new_reacs_file = open(fname+'.reac', 'w')
    for id_ in np.unique(kid_reduc):
        new_reacs_file.write(reacs_dict[id_])
    if add_He_reac:
        new_reacs_file.write('He + H2 = He + H2                                  0.00E+00     0.00        0    3 :                        : artificially added to involve He')
    new_reacs_file.close()
    if len(specs_lines) > 0:
        new_specs_file = open(fname+'.spec', 'w')

        try:
            if loss != -1 and imp_mols_flag:
                new_specs_file.write('! imp_loss: '+str(loss)+'\n')
                for i in range(1,3):
                    new_specs_file.write(specs_lines[i])
            else:
                for i in range(0,3):
                    new_specs_file.write(specs_lines[i])
            for m in np.unique(reduc_molecs):
                new_specs_file.write(specs_from_file[m])

            if not 'He' in np.unique(reduc_molecs):
                new_specs_file.write(' He           0.000E+00  -99999.000    -999.000   0.000000  0.000000    0     1     1    He  1')
        except Exception as e:
            print (e)
            print ('Error writing specs file:', specs_lines)

        
        new_specs_file.close()




