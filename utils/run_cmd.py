import subprocess

def run_cmd(cmd):
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result