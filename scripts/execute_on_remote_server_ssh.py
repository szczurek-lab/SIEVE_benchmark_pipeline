import argparse
import os
from typing import Union

import yaml
import paramiko


upload_files: dict = {
    "constants": "scripts/constants.py",
    "utils": "scripts/utils.py"
}


def parse_system_args() -> argparse.Namespace:
    sys_args_parser = argparse.ArgumentParser(description='Execute snakemake rules on a remote server with ssh.')
    sys_args_parser.add_argument('--config', help='path to configuration file', type=str)
    sys_args_parser.add_argument('--prefix', help='prefix relative to the root path on the remote server',
                                 metavar='FOLDER PREFIX', type=str)
    sys_args_parser.add_argument('--bash', help='bash script to run on the remote server',
                                 metavar='BASH SCRIPT', type=str)
    sys_args_parser.add_argument('--dataset', help='simulated dataset names', metavar='DATASET NAMES', type=str,
                                 nargs='+')
    sys_args_parser.add_argument('--upload', help='files to be uploaded (use relative path)',
                                 metavar='UPLOADED FILES', type=str, nargs='+')
    sys_args = sys_args_parser.parse_args()
    if not os.path.isfile(sys_args.config):
        raise ValueError('ERROR! File does not exist: ' + sys_args.config)
    upload_files["bash"] = sys_args.bash
    upload_files["upload"] = sys_args.upload
    return sys_args


def load_yaml(file_path: str) -> dict:
    with open(file_path, 'r') as fh:
        ssh_config: dict = yaml.load(fh, Loader=yaml.FullLoader)
    return ssh_config


def add_datasets_2_config(ssh_config: dict, datasets: list[str]):
    ssh_config['benchmark']['simulation']['datasetNames'] = datasets


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
    except paramiko.SSHException:
        print('SSH connection error.')
        exit(1)
    return jump_server_of_remote, remote_server


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
            upload_file_2_remote_server(remote_server,
                                        remote_server_sftp,
                                        val,
                                        os.path.join(root_path_on_server, val)
                                        )
        elif isinstance(val, list):
            for item in val:
                upload_file_2_remote_server(remote_server,
                                            remote_server_sftp,
                                            item,
                                            os.path.join(root_path_on_server, item)
                                            )


def execute_bash_script(remote_server: paramiko.SSHClient,
                        root_path_on_remote_server: str,
                        bash_file_name: str) -> None:
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


def close_servers(jump_server_of_remote: Union[None, paramiko.SSHClient], remote_server: paramiko.SSHClient):
    if jump_server_of_remote is not None:
        jump_server_of_remote.close()
    remote_server.close()


def main() -> None:
    sys_args: argparse.Namespace = parse_system_args()
    ssh_config = load_yaml(sys_args.config)
    add_datasets_2_config(ssh_config, sys_args.dataset)
    jump_server_of_remote, remote_server = create_ssh_connection(ssh_config)
    root_path_on_remote_server = os.path.join(ssh_config['servers']['remoteServer']['rootPath'], sys_args.prefix)
    upload_files_2_remote_server(remote_server,
                                 root_path_on_remote_server,
                                 sys_args.config.split('/')[-1],
                                 ssh_config
                                 )
    execute_bash_script(remote_server, root_path_on_remote_server, sys_args.bash)
    close_servers(jump_server_of_remote, remote_server)


if __name__ == '__main__':
    main()
