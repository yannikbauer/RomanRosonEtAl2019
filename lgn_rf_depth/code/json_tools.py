import json

def create_json_file(file_dir, data):
    # prepare json file for annotations
    with open(file_dir, 'w') as f:
        json.dump(data, f)


def read_from_json(file_name):
    with open(file_name) as f:
        data = json.load(f)
    return data