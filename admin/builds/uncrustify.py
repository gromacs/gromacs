import os
import subprocess

def do_build(context):
    os.environ['UNCRUSTIFY'] = context.env.get_uncrustify_command()
    warning_log = context.workspace.get_path_for_logfile('uncrustify.log', category='uncrustify')
    cmd = ['admin/uncrustify.sh', 'check', '--rev=HEAD^', '--warnings=' + warning_log]
    ret = subprocess.call(cmd)
    if ret == 1:
        with open(warning_log, 'r') as f:
            warnings = f.readlines()
        if len(warnings) <= 5:
            details = [x.rstrip() for x in warnings]
        else:
            uncrustify_count = 0
            cpyear_count = 0
            cpheader_count = 0
            for w in warnings:
                if 'needs uncrustify' in w:
                    uncrustify_count += 1
                if 'copyright year' in w:
                    cpyear_count += 1
                if 'copyright header' in w:
                    cpheader_count += 1
            details = []
            if uncrustify_count > 0:
                details.append('formatting issues in {0} files'.format(uncrustify_count))
            if cpyear_count > 0:
                details.append('copyright year missing in {0} files'.format(cpyear_count))
            if cpheader_count > 0:
                details.append('copyright header issues in {0} files'.format(cpheader_count))
        context.mark_unstable(reason='uncrustify.sh found issues', details=details)
    elif ret != 0:
        raise BuildError('uncrustify.sh failed to run')
