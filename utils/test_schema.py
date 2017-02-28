# Validates JSON input files for all benchmark tests against schema

import unittest
import os
import fnmatch
import json
import jsonschema

class JSONSchemaTestCase(unittest.TestCase):

    def __init__(self, filename, schema):
        super(JSONSchemaTestCase, self).__init__()
        self.filename = filename
        self.schema = schema

    def runTest(self):
        try:
            input = json.load(open(self.filename))
            jsonschema.validate(input, self.schema)
        except jsonschema.ValidationError:
            self.fail("JSON schema validation failed for file: " + self.filename)

def suite():
    suite = unittest.TestSuite()
    schema = json.load(open('input_schema.json'))
    path = "../test/benchmark/"
    filenames = [os.path.join(dirpath, f)
                 for dirpath, dirnames, files in os.walk(path)
                 for f in fnmatch.filter(files, '*.json')]
    suite.addTests(JSONSchemaTestCase(filename, schema) for filename in filenames)
    return suite

if __name__ == '__main__':

    unittest.TextTestRunner().run(suite())
