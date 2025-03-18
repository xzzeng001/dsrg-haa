import os

def clean_all_except_py():
    '''
    clean all files expect the python file for convenience
    '''

    os.system('rm -rf !(*.py|sub.sh)')
