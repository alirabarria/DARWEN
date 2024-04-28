import numpy as np
import pandas as pd
import subprocess
import os
import shutil
from pathlib import Path
from sklearn.metrics import r2_score, mean_absolute_percentage_error

imp_molecs = ['H2O', 'CH4', 'CO', 'CO2', 'NH3', 'HCN', 'C2H2', 'OH']
all_molecs = np.genfromtxt('1D_PAC/mols.txt', dtype=None, encoding=None)

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

def out_to_dfs(dir_, error_purposes=False, file_key='_out', header_keyword='!', return_only_mols=True):
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
            return None, errs
        else:
            return dfs, times, errs, irate, rvirate, rtype

    for file in matching_files:
        with open(os.path.join(dir_, file), 'r') as f:
            lines = f.readlines()
            # Find header line index based on the header keyword
            header_i = int(np.sum([header_keyword in line for i, line in enumerate(lines)])-1)

            # Extract error and conservation error
            err_line = lines[16].replace(' ', '')
            err = float(err_line.split('=')[1])
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
        return dfs[np.argmin(errs)], cons_err
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

def write_new_scheme(kid_reduc, reduc_molecs, fname='lira2023', verbose=False, loss=-1,
                    write_ids=False, imp_mols_flag=True, add_He_reac=False):
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
        new_reacs_file.write(reacs_from_file[id_])
    if add_He_reac:
        new_reacs_file.write('He + H2 = He + H2                                  0.00E+00     0.00        0    3 :                        : artificially added to involve He')
    new_reacs_file.close()
    new_specs_file = open(fname+'.spec', 'w')
    if loss != -1:
        if imp_mols_flag:
            new_specs_file.write('! imp_loss: '+str(loss)+'\n')
        else:
            new_specs_file.write('! loss: '+str(loss)+'\n')
        for i in range(1,3):
            new_specs_file.write(specs_lines[i])
    else:
        for i in range(0,3):
            new_specs_file.write(specs_lines[i])

    for m in np.unique(reduc_molecs):
        new_specs_file.write(specs_from_file[m])
    if not 'He' in np.unique(reduc_molecs):
        new_specs_file.write(' He           0.000E+00  -99999.000    -999.000   0.000000  0.000000    0     1     1    He  1')
    new_specs_file.close()



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

    # Run the bash command and capture the output
    result = subprocess.run(bash_command, shell=True, capture_output=True, text=True)

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

def extract_molecs(line, mode='array', link='->', check_line=0):
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
    sep = line.split('=')
    reacs = sep[0].strip().split('+')
    reacs = [r.strip() for r in reacs]
    prods = sep[1][:30].strip().split('+')
    prods = [p.strip() for p in prods]
    
    if mode == 'array' and check_line==0:
        return reacs, prods
    elif mode=='array' and check_line > 0:
        remainder = sep[1][30:30+check_line].strip().split('     ')
        return reacs, prods, remainder
    
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
def network_to_df(fname='1D_PAC/venot2020reduced.reac', start_i=11, mid_i=257, add_col='id', check_line=30):
    reacs_vreduc = open(fname, 'r')
    vreduc_lines = reacs_vreduc.readlines()
    i = start_i
    c = 0
    columns = ['reacs', 'prods', 'molecs']
    if check_line > 0:
        columns=np.append(columns, ['to_check'])
        
    if add_col != '':
        reacs_df = pd.DataFrame(columns=np.append(columns, [add_col]))
    else:
        reacs_df = pd.DataFrame(columns=columns)
        
    while i < len(vreduc_lines):
        if check_line==0:
            reacs, prods = extract_molecs(vreduc_lines[i], check_line=check_line)
        else:
            reacs, prods, remainder = extract_molecs(vreduc_lines[i], check_line=check_line)
        reacs_df.at[c, 'reacs'] = reacs
        reacs_df.at[c, 'prods'] = prods
        reacs_df.at[c, 'molecs'] = np.append(reacs, prods)
        if check_line > 0:
            reacs_df.at[c, 'to_check'] = remainder
        if i < mid_i:
            c+=1
            i+=1

        else:
            i+=1
            if i < len(vreduc_lines):
                while len(vreduc_lines[i][:4].replace(' ', '')) == 0:

                    i+=1      
                c+=1
    reacs_df.index = range(len(reacs_df))
    return reacs_df


reacs_in_venot = network_to_df()
reacs_in_full = network_to_df(fname='1D_PAC/venot2020.reac', start_i=0, mid_i=837, add_col='')
reacs_file = open('1D_PAC/venot2020.reac', 'r')
reacs_lines = reacs_file.readlines()
reacs_from_file = {}
i = 0
c = 0
while i < len(reacs_lines):
    if i < 837:
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
    reacs, prods, remainder = extract_molecs(reacs_from_file[c-1], check_line=30)
    for ind, row in reacs_in_venot.iterrows():
        if np.all(reacs == row['reacs']) and np.all(prods == row['prods']):
            if np.all(remainder == row['to_check']):
                reacs_in_venot.at[ind, 'id'] = c - 1
            else:
                print ('signaling', ind, c - 1, remainder)
#reacs_in_venot = reacs_in_venot.dropna()

specs_file = open('1D_PAC/venot2020.spec', 'r')
specs_lines = specs_file.readlines()
specs_from_file = {}
for i in range(3,len(specs_lines)):
    sp = specs_lines[i].split(' ')[1]
    specs_from_file[sp] = specs_lines[i]

for i in range(len(reacs_from_file)):
    arr = reacs_from_file[i].split(':')
    splitted = arr[0].split()
    rtype = splitted[-1]
    gamma = splitted[-2]
    beta = splitted[-3]
    alpha = splitted[-4]
    reacs_in_full.at[i, 'rtype'] = int(rtype)
    reacs_in_full.at[i, 'alpha'] = np.float64(alpha)
    reacs_in_full.at[i, 'beta'] = np.float64(beta)
    reacs_in_full.at[i, 'gamma'] = np.float64(gamma)

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

    with open(output_file, 'w') as output_file:
        process = subprocess.Popen(command, shell=True, stdout=output_file, stderr=subprocess.STDOUT, 
                                       cwd=out_path)
        process.wait()  # Wait for the process to finish


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


