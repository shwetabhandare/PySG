import os
import generate_model
from generate_model import *


testFile = sys.argv[1]
a_count, t_count, c_count, g_count, a_repeated_count, t_repeated_count, au_count = generate_model.get_ATGC_counts(testFile);
print a_count
normalized_a_count = generate_model.get_normalized_array(a_count);
