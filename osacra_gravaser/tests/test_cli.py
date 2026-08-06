import subprocess
import sys
import os

def run_cli_test(mode, extra_args=[]):
    print(f"\n--- Testing CLI Mode: {mode.upper()} ---")
    cli_script = os.path.join(os.path.dirname(__file__), '..', 'cli.py')
    cmd = [sys.executable, cli_script, "--mode", mode, "--render", "none", "--steps", "50"]
    cmd.extend(extra_args)
    
    print(f"Running command: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"SUCCESS! Output:\n{result.stdout.strip()[-250:]}") # Print the last 250 chars of output to verify it finished
    except subprocess.CalledProcessError as e:
        print(f"FAILED! Error output:\n{e.stderr}")
        sys.exit(1)

def main():
    print("=========================================================")
    print(" CEFE CLI Automated Headless Tests")
    print("=========================================================")
    
    run_cli_test("qm")
    run_cli_test("holographic")
    run_cli_test("qft")
    
    print("\nAll CLI tests passed successfully!")

if __name__ == "__main__":
    main()
