import subprocess
import os


# NOTE: Python 3.6. Commented for 3.7
def get_git_commit():
    """Get the latest commit hash"""
    completed_process = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                                       #capture_output=True, text=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       cwd=os.path.dirname(os.path.abspath(__file__)))
    git_commit = completed_process.stdout.decode('ascii').strip()
    return git_commit


# NOTE: Python 3.6. Commented for 3.7
def get_git_tag():
    """Get the current git tag for the repository."""
    completed_process = subprocess.run(['git', 'describe', '--tags'],
                                       #capture_output=True, text=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       cwd=os.path.dirname(os.path.abspath(__file__)))
    git_tag = completed_process.stdout.decode('ascii').strip()  # Strip probably not needed
    return git_tag

# NOTE: Python 3.6. Commented for 3.7
# NOTE: This only checks for modified tracked files, not untracked files.
def is_branch_clean():
    """Return whether current branch/commit is dirty"""
    # NOTE: Flag --always enables fallback on latest commit if no tag available
    completed_process = subprocess.run(['git', 'describe', '--tags', '--always', '--dirty'],
                                       #capture_output=True, text=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       cwd=os.path.dirname(os.path.abspath(__file__)))
    git_tag_or_hash = completed_process.stdout.decode('ascii').strip()
    return 'dirty' not in git_tag_or_hash

# NOTE: Python 3.6. Commented for 3.7
def get_git_reponame():
    """Get the name of the git repository"""
    completed_process = subprocess.run(['git', 'config', '--get', 'remote.origin.url'],
                                       #capture_output=True, text=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       cwd=os.path.dirname(os.path.abspath(__file__)))
    repo_name = completed_process.stdout.decode('ascii').strip()  # Strip probably not needed

    #Remove trailing and leading information
    repo_name = os.path.basename(repo_name)[:-4]

    return repo_name
