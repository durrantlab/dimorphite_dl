# Copyright 2018 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Main entry script for running protonation.
"""

import argparse

import protonation_functions as pf

parser = argparse.ArgumentParser(description="Protonates small moleucles.")
parser.add_argument('--min_ph', metavar='MIN', type=float, default=6.4,
                    help='Minimum pH to consider.')
parser.add_argument('--max_ph', metavar='MAX', type=float, default=8.4,
                    help='Maximum pH to consider.')
parser.add_argument('--st_dev', metavar='STD', type=float, default=1.0,
                    help='Standard devation range (number of standard devs).')
parser.add_argument('--smiles', type=str,
                    help='SMILE string to protonate.')
parser.add_argument('--smiles_file', type=str,
                    help='File which contains SMILES strings to protonate.')
parser.add_argument('--output_file', type=str,
                    help='File to write protonated SMILES. (Optional)')


if __name__ == "__main__":
    args = vars(parser.parse_args())

    output = pf.protonate(args)

    if 'output_file' in args:
        with open(args['output_file'], 'w') as file:
            file.write("\n".join(output))
    else:
        for out in output:
            print(out)
