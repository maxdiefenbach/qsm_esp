import numpy as np
import SimpleITK as sitk
from utils.utils import *
import os

os.environ['SITK_SHOW_COMMAND'] = \
    '/Applications/ITK-SNAP.app/Contents/MacOS/ITK-SNAP'

close_all()


f_cosmos = '../data/20170327_qsm2016_recon_challenge/data/chi_cosmos.nii.gz'
f_chi33 = '../data/20170327_qsm2016_recon_challenge/data/chi_33.nii.gz'

cosmos = sitk.ReadImage((f_cosmos))
chi33 = sitk.ReadImage((f_chi33))

cosmos.SetDirection(np.abs(cosmos.GetDirection()))
chi33.SetDirection(np.abs(chi33.GetDirection()))

sitk.Show(cosmos, 'cosmos')
sitk.Show(chi33, 'chi33')

# cosmos_flipped = sitk.Flip(cosmos, [False, True, False], True)
# chi33_flipped = sitk.Flip(chi33, [False, True, False], True)

# sitk.Show(cosmos_flipped, 'cosmos_flipped')
# sitk.Show(chi33_flipped, 'chi33_flipped')

sitk.WriteImage(cosmos, 'cosmos_flipped.nii.gz')
sitk.WriteImage(chi33, 'chi33_flipped.nii.gz')
