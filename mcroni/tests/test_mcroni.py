import pytest
import subprocess
import os

import mcroni.seqFunctions as sf

# Data directory
#data_dir = dir+'/../data'

def test_classify_variant():
      with open('fake.fa', 'w') as f: # make a fake fasta file
          f.write('>FAKE_1\nAAA\nFAKE_2\nTTT')
      assert sf.classify_variant('AAA', variants_db='fake.fa')=='FAKE_1'
