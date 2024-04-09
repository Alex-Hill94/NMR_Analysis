import os
import csv

def process_acqus_file(folder, filename, linestart='##$NS'):
    file_path = os.path.join(folder, filename)
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith(linestart) and line[len(linestart)] in ('', ' ', '='):
                    _, value = line.split('=')
                    return value.strip()
        return None
    except Exception as e:
        print(f"Error processing file '{filename}': {e}")
        return None



def extract_root_and_experiment(folder):
    # Split folder path into Root and Experiment
    root, experiment = os.path.split(folder)
    return root, experiment
    
def search_and_save_to_csv(root_folder='.', output_csv='output.csv'):
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Root', 'Sample', 'Experiment', 'Filename', 'PULPROG', 'NS (Num. Scans)', 'TD0 (Loop count)', 'Total Scans', 'RG', 'TD (Size of FID)', 'Time (seconds)']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for foldername, subfolders, filenames in os.walk(root_folder):
            for filename in filenames:
                if filename.lower() == 'acqus':
                    ns_value = process_acqus_file(foldername, filename, linestart='##$NS')
                    td0_value = process_acqus_file(foldername, filename, linestart='##$TD0')
                    td_value = process_acqus_file(foldername, filename, linestart='##$TD')
                    pulprog_value = process_acqus_file(foldername, filename, linestart='##$PULPROG')
                    date_value = process_acqus_file(foldername, filename, linestart='##$DATE')
                    rg_value = process_acqus_file(foldername, filename, linestart='##$RG')
                    date_st_value = process_acqus_file(foldername, filename, linestart='##$DATE_START')

                    if ns_value is not None:
                        root_pr, experiment = extract_root_and_experiment(foldername)
                        root, sample = extract_root_and_experiment(root_pr)
                        writer.writerow({
                            'Root': root,
                            'Sample': sample,
                            'Experiment': experiment,
                            'Filename': filename,
                            'PULPROG': pulprog_value,
                            'NS (Num. Scans)': ns_value,
                            'TD0 (Loop count)': td0_value,
                            'Total Scans': int(ns_value) * int(td0_value),
                            'RG': rg_value,
                            'TD (Size of FID)': td_value,
                            'Time (seconds)': int(date_value) - int(date_st_value)
                        })

# Replace '.' with the root folder path if you want to start the search from a specific directory
# Replace 'output.csv' with the desired CSV file name
search_and_save_to_csv('/Users/alexhill/Desktop/Metabolomics/Rerun_Data/20240119_RERUN_3', '20240119_RERUN_3.csv')