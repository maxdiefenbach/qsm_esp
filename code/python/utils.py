import os
import psutil
import numpy as np
import SimpleITK as sitk


def close_all(viewer=['JavaApplicationStub', 'ITK-SNAP']):
    for v in viewer:
        kill_processes(v)


def killITK():
    kill_processes('ITK-SNAP')


def killIJ():
    kill_processes('JavaApplicationStub')


def kill_processes(procname):
    for proc in psutil.process_iter():
        pinfo = proc.as_dict(attrs=['pid', 'name'])
        if pinfo['name'] is None:
            continue
        if procname in pinfo['name']:
            proc.kill()
            print('killed:', pinfo)


def show_arr3d(arr, tool='itk-snap', title=None):
    SITK_SHOW_COMMAND = os.environ.get('SITK_SHOW_COMMAND')

    print(tool)
    if tool.lower() == 'itk-snap':
        os.environ['SITK_SHOW_COMMAND'] = \
            '/Applications/ITK-SNAP.app/Contents/MacOS/ITK-SNAP'
    elif tool.lower() == 'imagej':
        os.environ['SITK_SHOW_COMMAND'] = \
            '/Applications/ImageJ/ImageJ.app/Contents/MacOS/JavaApplicationStub'
    elif tool.lower() == 'slicer':
        os.environ['SITK_SHOW_COMMAND'] = \
            '/Applications/Slicer.app/Contents/MacOS/Slicer'
    else:
        print('tool not known.')
        os.environ['SITK_SHOW_COMMAND'] = tool

    # SimpleITK takes array's in scikit.ndimage convention [z, y, x] and
    # ITK-snap flips i and j directions
    # therefore transpose array to [k, j, i] and flip
    # and then flip i, j to [k, i, j] => transpose order [2, 0, 1]
    if title is not None:
        sitk.Show(sitk.GetImageFromArray(
            np.transpose(arr.astype(float), [2, 0, 1])),
                  title)
    else:
        sitk.Show(sitk.GetImageFromArray(
            np.transpose(arr.astype(float), [2, 0, 1])))


    if SITK_SHOW_COMMAND is not None:
        os.environ['SITK_SHOW_COMMAND'] = SITK_SHOW_COMMAND


def read_array(filename):
    return np.transpose(sitk.GetArraryFromImage(
        sitk.ReadImage(filename)),
                        [2, 0, 1])


def write_array(array, filename):
    return sitk.WriteImage(sitk.GetImageFromArray(
        np.transpose(array, [2, 0, 1])), filename)
