import subprocess

def run_cmd(cmd, input_data=None, capture_output=True, text=True, check=True, shell=False, **kwargs):
    result = subprocess.run(cmd, 
                            input=input_data, 
                            capture_output=capture_output, 
                            text=text, 
                            check=check,
                            shell=shell,
                            **kwargs)
    if result.returncode != 0:
        print(f"Command failed with return code {result.returncode}: {cmd}")
        print(f"stderr: {result.stderr}")
        raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)
    return result