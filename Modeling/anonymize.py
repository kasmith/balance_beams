import os
import pandas as pd
import json

# Ensure that directories exist
os.makedirs('anonymized_output', exist_ok=True)
os.makedirs('anonymized_output/raw_data', exist_ok=True)

# A dictionary to map WIDs to their anonymized versions
WID_to_anonymized = {}

def anonymize_wid(wid):
    """Anonymize a WID, using the WID_to_anonymized dictionary to ensure consistency across files."""
    if wid not in WID_to_anonymized:
        WID_to_anonymized[wid] = f"S{len(WID_to_anonymized)+1:03}"
    return WID_to_anonymized[wid]

# List of CSV files to process
files = [
    'models/geomat/BB_GeoMatData.csv',
    'models/ferretti/BB_FerrettiData.csv',
    'models/combined/BB_CombData.csv',
    'models/learn_benefit/BB_BeneData.csv'
]

# Process each CSV file
for file in files:
    # Load the file into a pandas DataFrame
    df = pd.read_csv(file)

    # If the DataFrame has a 'WID' column, anonymize it
    if 'WID' in df.columns:
        df['WID'] = df['WID'].map(anonymize_wid)

    # Determine the path to write the anonymized file to
    filename = os.path.basename(file)
    new_filepath = os.path.join('anonymized_output/raw_data', filename)

    # Write the anonymized DataFrame to a new CSV file
    df.to_csv(new_filepath, index=False)

# Process the specific JSON file
json_file_path = 'output/learn_bene/base_strats_joint_percept_params.json'
with open(json_file_path, 'r') as json_file:
    data = json.load(json_file)

# If the JSON object has an 'ind_strategies' key, anonymize its subkeys
if 'ind_strategies' in data:
    data['ind_strategies'] = {anonymize_wid(k): v for k, v in data['ind_strategies'].items()}# List of JSON files to process
json_files = [
    'output/learn_bene/base_strats_joint_percept_params.json',
    'output/comb_strats/base_strats_all_params.json'
]

# Process each JSON file
for file in json_files:
    with open(file, 'r') as json_file:
        data = json.load(json_file)

    # If the JSON object has an 'ind_strategies' key, anonymize its subkeys
    if 'ind_strategies' in data:
        data['ind_strategies'] = {anonymize_wid(k): v for k, v in data['ind_strategies'].items()}

    # Determine the path to write the anonymized file to
    new_dirpath = file.replace('output', 'anonymized_output', 1)
    os.makedirs(os.path.dirname(new_dirpath), exist_ok=True)

    # Write the anonymized JSON data to a new file
    with open(new_dirpath, 'w') as json_file:
        json.dump(data, json_file)


# Continue with the directory walk as before
for dirpath, dirnames, filenames in os.walk('output'):
    for filename in filenames:
        if filename.endswith('.csv'):
            # Load the file into a pandas DataFrame
            df = pd.read_csv(os.path.join(dirpath, filename))

            # If the DataFrame has a 'WID' column, anonymize it
            if 'WID' in df.columns:
                df['WID'] = df['WID'].map(anonymize_wid)

            # Determine the path to write the anonymized file to
            new_dirpath = dirpath.replace('output', 'anonymized_output', 1)
            os.makedirs(new_dirpath, exist_ok=True)
            new_filepath = os.path.join(new_dirpath, filename)

            # Write the anonymized DataFrame to a new CSV file
            df.to_csv(new_filepath, index=False)
