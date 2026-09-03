import yaml
import json
with   open("reims_schema.yaml", 'r') as yaml_in, \
       open("reims_schema.json", "w") as json_out:
    yaml_object = yaml.safe_load(yaml_in)
    json.dump(yaml_object, json_out, indent=2)