import argparse
import os
import time
from typing import Union

import yaml
import paramiko
from paramiko.agent import AgentRequestHandler

upload_files: dict[str, str] = {
    "constants": "scripts/constants.py",
    "utils": "scripts/utils.py"
}

WAITING_TIME_1: int = 10
WAITING_TIME_2: int = 60


def parse_system_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(description='Execute snakemake rules on a remote server with git.')
    sys_args_parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='relative path to configuration file'
    )
    sys_args_parser.add_argument(
        '--bash',
        metavar='BASH SCRIPT',
        type=str,
        required=True,
        help='relative path to bash script to run on the remote server'
    )
    sys_args_parser.add_argument(
        '--upload',
        metavar='UPLOADED FILES',
        type=str,
        nargs='+',
        required=True,
        help='relative path to files to be uploaded'
    )
    sys_args_parser.add_argument(
        '--dataset',
        metavar='DATASET NAMES',
        type=str,
        nargs='+',
        required=True,
        help='simulated dataset names'
    )
    sys_args_parser.add_argument(
        '--datatype',
        metavar='DATA TYPE NAMES',
        type=str,
        nargs='+',
        required=True,
        help='simulated data type names'
    )
    sys_args_parser.add_argument(
        '--snkmkdir',
        metavar='PATH TO SNAKEMAKE',
        type=str,
        required=True,
        help='the folder where snakemake rules will be executed, relative to the root path of the remote server'
    )
    sys_args_parser.add_argument(
        '--gitdir',
        metavar='GIT REPO CLONED',
        type=str,
        required=True,
        help='the folder where git repository is cloned, relative to the root path of the remote server'
    )
    sys_args_parser.add_argument(
        '--repopath',
        metavar='GIT REPO',
        type=str,
        required=True,
        help='the absolute path to the git repo'
    )
    sys_args = sys_args_parser.parse_args()

    if not os.path.isfile(sys_args.config):
        raise ValueError('ERROR! File does not exist: ' + sys_args.config)

    if not os.path.isfile(sys_args.bash):
        raise ValueError('ERROR! File does not exist: ' + sys_args.bash)

    upload_files["bash"] = sys_args.bash
    upload_files["upload"] = sys_args.upload
    return sys_args


def load_yaml(file_path: str) -> dict:
    with open(file_path, 'r') as fh:
        ssh_config: dict = yaml.load(fh, Loader=yaml.FullLoader)
    return ssh_config


def add_data_2_config(ssh_config: dict, datasets: list[str], datatypes: list[str]):
    ssh_config['benchmark']['simulation']['datasetNames'] = datasets
    ssh_config['benchmark']['simulation']['dataTypeNames'] = datatypes


def create_ssh_connection(ssh_config: dict) -> \
        tuple[Union[None, paramiko.SSHClient], paramiko.SSHClient]:
    jump_server_of_remote = None

    remote_server: paramiko.SSHClient = paramiko.SSHClient()
    remote_server.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    try:
        if ssh_config['servers']['jumpServerOfRemote']['required']:
            jump_server_of_remote = paramiko.SSHClient()
            jump_server_of_remote.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            jump_server_of_remote.connect(
                hostname=ssh_config['servers']['jumpServerOfRemote']['host'],
                port=ssh_config['servers']['jumpServerOfRemote']['port'],
                username=ssh_config['servers']['jumpServerOfRemote']['user'],
                password=ssh_config['servers']['jumpServerOfRemote']['password']
            )
            j2r_channel = jump_server_of_remote.get_transport().open_channel(
                kind='direct-tcpip',
                dest_addr=(ssh_config['servers']['remoteServer']['host'],
                           ssh_config['servers']['remoteServer']['port']),
                src_addr=(ssh_config['servers']['jumpServerOfRemote']['host'],
                          ssh_config['servers']['jumpServerOfRemote']['port'])
            )
            remote_server.connect(
                hostname=ssh_config['servers']['remoteServer']['host'],
                port=ssh_config['servers']['remoteServer']['port'],
                username=ssh_config['servers']['remoteServer']['user'],
                password=ssh_config['servers']['remoteServer']['password'],
                sock=j2r_channel
            )
        else:
            remote_server.connect(
                hostname=ssh_config['servers']['remoteServer']['host'],
                port=ssh_config['servers']['remoteServer']['port'],
                username=ssh_config['servers']['remoteServer']['user'],
                password=ssh_config['servers']['remoteServer']['password']
            )
    except paramiko.SSHException as e:
        print(e)
        print('SSH connection error.')
        exit(1)

    print('Remote server connected.')
    return jump_server_of_remote, remote_server


def get_full_git_repo(ssh_config: dict, repo_path: str) -> str:
    full_path: str = 'ssh://' + ssh_config['servers']['localServer']['user'] + '@' + \
                     ssh_config['servers']['localServer']['host'] + \
                     ':' + str(ssh_config['servers']['localServer']['port'])
    full_path += os.path.abspath(repo_path)
    return full_path


def clone_git_repo(
        ssh_config: dict,
        ssh_server: paramiko.SSHClient,
        sftp_server: paramiko.SFTPClient,
        remote_git_dir: str,
        local_git_repo: str
) -> None:
    # get repo name
    comp: list[str] = local_git_repo.strip().split('/')
    repo_name: str = comp[-1] if comp[-1].strip() != '' else comp[-2]

    # use an interactive shell to emulate a login session
    channel = ssh_server.get_transport().open_channel(kind='session')
    channel.get_pty()
    channel.invoke_shell()
    AgentRequestHandler(channel)

    try:
        print(sftp_server.stat(os.path.join(remote_git_dir, repo_name)))
        print('Git repo detected. Performing \'git pull\'...')

        channel.send('cd ' + os.path.join(remote_git_dir, repo_name) + '\n')
        if channel.recv_ready():
            print(channel.recv(1024).decode('utf-8'))

        channel.send('git pull origin ' + ssh_config['benchmark']['git']['branchName'] + '\n')
        time.sleep(WAITING_TIME_1)
        if channel.recv_ready():
            print(channel.recv(1024).decode('utf-8'))

        channel.send(ssh_config['servers']['localServer']['password'] + '\n')
        time.sleep(WAITING_TIME_2)
        if channel.recv_ready():
            print(channel.recv(9999).decode('utf-8'))

        print('Git repo updated.')
    except IOError:
        print('Git repo does not exist. Performing \'git clone\'...')

        channel.send('mkdir -p ' + remote_git_dir + '\n')
        if channel.recv_ready():
            print(channel.recv(1024).decode('utf-8'))

        channel.send('cd ' + remote_git_dir + '\n')
        if channel.recv_ready():
            print(channel.recv(1024).decode('utf-8'))

        channel.send('git clone --single-branch --branch ' + ssh_config['benchmark']['git']['branchName'] + ' ' + local_git_repo + '\n')
        time.sleep(WAITING_TIME_1)
        if channel.recv_ready():
            print(channel.recv(1024).decode('utf-8'))

        channel.send(ssh_config['servers']['localServer']['password'] + '\n')
        time.sleep(WAITING_TIME_2)
        if channel.recv_ready():
            print(channel.recv(9999).decode('utf-8'))

        print('Git repo clone finished.')


def create_sftp(server: paramiko.SSHClient) -> paramiko.SFTPClient:
    return server.open_sftp()


def dump_config_on_remote_server(sftp: paramiko.SFTPClient, file_name: str, ssh_config: dict) -> None:
    try:
        print(sftp.stat(file_name))
        print(file_name + ' exists on remote server. Skip creation.')
        return
    except IOError:
        print('Creating ' + file_name + '...')
        with sftp.open(file_name, 'w') as fh:
            yaml.dump(ssh_config, fh)


def upload_file_2_remote_server(ssh: paramiko.SSHClient, sftp: paramiko.SFTPClient,
                                local_path: str, remote_path: str) -> None:
    try:
        print(sftp.stat(remote_path))
        print(remote_path + ' exists on remote server. Skip copy.')
        return
    except IOError:
        print('Copying to ' + remote_path + '...')
        ssh.exec_command('mkdir -p ' + os.path.dirname(remote_path))
        print(local_path)
        sftp.put(local_path, remote_path)


def upload_files_2_remote_server(
        remote_server: paramiko.SSHClient,
        root_path_on_server: str,
        config_file_name: str,
        ssh_config: dict
) -> None:
    remote_server_sftp: paramiko.SFTPClient = create_sftp(remote_server)
    remote_server.exec_command("mkdir -p " + os.path.dirname(root_path_on_server))
    dump_config_on_remote_server(remote_server_sftp, os.path.join(root_path_on_server, config_file_name), ssh_config)
    for key, val in upload_files.items():
        if isinstance(val, str):
            upload_file_2_remote_server(
                remote_server,
                remote_server_sftp,
                val,
                os.path.join(root_path_on_server, val)
            )
        elif isinstance(val, list):
            for item in val:
                upload_file_2_remote_server(
                    remote_server,
                    remote_server_sftp,
                    item,
                    os.path.join(root_path_on_server, item)
                )


def execute_bash_script(
        remote_server: paramiko.SSHClient,
        root_path_on_remote_server: str,
        bash_file_name: str
) -> None:
    print('Changing working directory on remote server to ' + root_path_on_remote_server +
          ' and executing ' + bash_file_name)

    # use an interactive shell to emulate a login session
    channel = remote_server.get_transport().open_channel(kind='session')
    channel.get_pty()
    channel.invoke_shell()
    channel.send('cd ' + root_path_on_remote_server + '\n')
    channel.send('source ' + bash_file_name + '\n')
    exit_status = channel.recv_exit_status()
    if exit_status == 0:
        print('Succeed!')
    else:
        print('Something is wrong. Please check the following file on remote server: ' + root_path_on_remote_server +
              '/out')


def push_commits(
        remote_server: paramiko.SSHClient,
        ssh_config: dict,
        remote_git_dir: str,
        local_git_repo: str
) -> None:
    # get repo name
    comp: list[str] = local_git_repo.strip().split('/')
    repo_name: str = comp[-1] if comp[-1].strip() != '' else comp[-2]

    channel = remote_server.get_transport().open_channel(kind='session')
    channel.get_pty()
    channel.invoke_shell()

    channel.send('cd ' + os.path.join(remote_git_dir, repo_name) + '\n')
    if channel.recv_ready():
        print(channel.recv(1024).decode('utf-8'))

    channel.send('git push -u origin ' + ssh_config['benchmark']['git']['branchName'] + '\n')
    time.sleep(WAITING_TIME_1)
    if channel.recv_ready():
        print(channel.recv(1024).decode('utf-8'))

    channel.send(ssh_config['servers']['localServer']['password'] + '\n')
    time.sleep(WAITING_TIME_2)
    if channel.recv_ready():
        print(channel.recv(9999).decode('utf-8'))


def close_servers(jump_server_of_remote: Union[None, paramiko.SSHClient], remote_server: paramiko.SSHClient):
    if jump_server_of_remote is not None:
        jump_server_of_remote.close()
    remote_server.close()


def main() -> None:
    sys_args: argparse.Namespace = parse_system_args()

    # load and update configuration
    ssh_config = load_yaml(sys_args.config)
    add_data_2_config(ssh_config, sys_args.dataset, sys_args.datatype)

    # establish connection
    jump_server_of_remote, remote_server = create_ssh_connection(ssh_config)

    remote_git_dir = os.path.join(
            ssh_config['servers']['remoteServer']['rootPath'],
            sys_args.gitdir
    )
    local_git_repo = get_full_git_repo(
            ssh_config,
            sys_args.repopath
    )

    # clone the git repository on remote server
    clone_git_repo(
        ssh_config,
        remote_server,
        create_sftp(remote_server),
        remote_git_dir,
        local_git_repo
    )

    # upload snakemake-related files
    root_path_on_remote_server = os.path.join(ssh_config['servers']['remoteServer']['rootPath'], sys_args.snkmkdir)
    upload_files_2_remote_server(
        remote_server,
        root_path_on_remote_server,
        sys_args.config.split('/')[-1],
        ssh_config
    )
    execute_bash_script(
        remote_server,
        root_path_on_remote_server,
        sys_args.bash
    )

    # push committed file changes
    push_commits(
        remote_server,
        ssh_config,
        remote_git_dir,
        local_git_repo
    )

    # close connection
    close_servers(jump_server_of_remote, remote_server)


if __name__ == '__main__':
    main()
