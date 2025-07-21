import sys
import os
import fcsparser

def get_event_count(fcs_file):
    try:
        meta = fcsparser.parse(fcs_file, meta_data_only=True, reformat_meta=False)
        return int(meta['$TOT'])
    except TypeError:
        try:
            meta, _ = fcsparser.parse(fcs_file, reformat_meta=False)
            return int(meta['$TOT'])
        except Exception as e:
            print(f"Error reading FCS file {fcs_file}: {e}")
            return 0
    except Exception as e:
        print(f"Error reading FCS file {fcs_file}: {e}")
        return 0

def process_directory(directory):
    total_events = 0
    fcs_files_processed = 0

    for root, _, files in os.walk(directory):
        for file in files:
            if file.lower().endswith('.fcs'):
                fcs_path = os.path.join(root, file)
                events = get_event_count(fcs_path)
                total_events += events
                fcs_files_processed += 1

    return total_events, fcs_files_processed

def main():
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <path_to_directory>")
        sys.exit(1)

    directory = sys.argv[1]
    
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory")
        sys.exit(1)

    total_events, files_processed = process_directory(directory)

    print(f"\nTotal number of events across all FCS files: {total_events}")
    print(f"Number of FCS files processed: {files_processed}")

if __name__ == "__main__":
    main()