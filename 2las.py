import os
import laspy
import numpy as np
from plyfile import PlyData, PlyElement

def load_ply2las(path):
    plydata = PlyData.read(path)

    xyz = np.stack((np.asarray(plydata.elements[0]["x"]),
                    np.asarray(plydata.elements[0]["y"]),
                    np.asarray(plydata.elements[0]["z"])),  axis=1)
    opacities = np.asarray(plydata.elements[0]["opacity"])[..., np.newaxis]

    features_dc = np.zeros((xyz.shape[0], 3, 1))
    features_dc[:, 0, 0] = np.asarray(plydata.elements[0]["f_dc_0"])
    features_dc[:, 1, 0] = np.asarray(plydata.elements[0]["f_dc_1"])
    features_dc[:, 2, 0] = np.asarray(plydata.elements[0]["f_dc_2"])

    # extra_f_names = [p.name for p in plydata.elements[0].properties if p.name.startswith("f_rest_")]
    # extra_f_names = sorted(extra_f_names, key = lambda x: int(x.split('_')[-1]))
    # assert len(extra_f_names)==3*(self.max_sh_degree + 1) ** 2 - 3
    # features_extra = np.zeros((xyz.shape[0], len(extra_f_names)))
    # for idx, attr_name in enumerate(extra_f_names):
    #     features_extra[:, idx] = np.asarray(plydata.elements[0][attr_name])
    # # Reshape (P,F*SH_coeffs) to (P, F, SH_coeffs except DC)
    # features_extra = features_extra.reshape((features_extra.shape[0], 3, (self.max_sh_degree + 1) ** 2 - 1))

    scale_names = [p.name for p in plydata.elements[0].properties if p.name.startswith("scale_")]
    scale_names = sorted(scale_names, key = lambda x: int(x.split('_')[-1]))
    scales = np.zeros((xyz.shape[0], len(scale_names)))
    for idx, attr_name in enumerate(scale_names):
        scales[:, idx] = np.asarray(plydata.elements[0][attr_name])

    rot_names = [p.name for p in plydata.elements[0].properties if p.name.startswith("rot")]
    rot_names = sorted(rot_names, key = lambda x: int(x.split('_')[-1]))
    rots = np.zeros((xyz.shape[0], len(rot_names)))
    for idx, attr_name in enumerate(rot_names):
        rots[:, idx] = np.asarray(plydata.elements[0][attr_name])


     # 1. Create a Las
    header = laspy.LasHeader(point_format=2, version="1.4")

    for i in range(features_dc.shape[1]):
        header.add_extra_dim(laspy.ExtraBytesParams(name='f_dc_{}'.format(i), type=np.float32))

    # for i in range(f_rest.shape[1]):
    #     header.add_extra_dim(laspy.ExtraBytesParams(name='f_rest_{}'.format(i), type=np.float32))

    header.add_extra_dim(laspy.ExtraBytesParams(name='opacity', type=np.float32))

    for i in range(scales.shape[1]):
        header.add_extra_dim(laspy.ExtraBytesParams(name='scale_{}'.format(i), type=np.float32))

    for i in range(rots.shape[1]):
        header.add_extra_dim(laspy.ExtraBytesParams(name='rot_{}'.format(i), type=np.float32))

    # 2. Create a Las
    out_las = laspy.LasData(header)

    # 3. Fill the Las
    out_las.x = xyz[:, 0]
    out_las.y = xyz[:, 1]
    out_las.z = xyz[:, 2]
    # Fill the extra attributes
    for i in  range(features_dc.shape[1]):
        setattr(out_las, f"f_dc_{i}", features_dc[:, i, 0])
    
    # for i in  range(f_rest.shape[1]):
    #     setattr(out_las, f"f_rest_{i}", f_rest[:, i])
    
    out_las.opacity = opacities[:, 0]

    for i in  range(scales.shape[1]):
        setattr(out_las, f"scale_{i}", scales[:, i])

    for i in  range(rots.shape[1]):
        setattr(out_las, f"rot_{i}", rots[:, i])

    # Save the LAS file
    out_las.write(path.replace(".ply", ".las"))


if __name__ == "__main__":
    load_ply2las("./model/point_cloud.ply")