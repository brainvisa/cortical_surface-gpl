from __future__ import absolute_import
from brainvisa.processes import *
from six.moves import range
try:
    from brainvisa import anatomist
except ImportError:
    pass

name = 'Flat-map image from rectangular mesh'

def validation():
    try:
        from brainvisa import anatomist
    except ImportError:
        raise ValidationError('Anatomist is needed')

signature = Signature(
    'flat_texture', ReadDiskItem('Texture', 'aims texture formats'),
    'flat_mesh', ReadDiskItem('Rectangular flat mesh',
                              'anatomist mesh formats'),
    'image_2d', WriteDiskItem('2D Image', 'aims writable volume formats'),
    'width', Integer(),
    'height', Integer(),
    'rgb_interpolation', Boolean(),
    'output_rgb_image', Boolean(),
    'palette', String(),
    'background_value', Float(),
    'keep_aspect_ratio', Boolean(),
)


def initialization(self):
    self.linkParameters('flat_mesh', 'flat_texture')
    self.width = 512
    self.height = 0
    self.rgb_interpolation = True
    self.output_rgb_image = False
    self.palette = 'B-W LINEAR'
    self.background_value = 0
    self.keep_aspect_ratio = True


def execution_mainthread(self, context):

    from soma.qt_gui import qt_backend
    import numpy as np

    a = anatomist.Anatomist()
    amesh = a.loadObject(self.flat_mesh)
    amesh.setMaterial(diffuse=[1., 1., 1., 1.], lighting=False)
    tex = aims.read(self.flat_texture.fullPath())
    atex = a.toAObject(tex)
    atex.releaseAppRef()
    if self.rgb_interpolation:
        interp = 'rgb'
    else:
        interp = 'palette'
    a.execute('TexturingParams', objects=[atex], interpolation=interp)
    ntex = np.asarray(tex[0])
    tmin = np.min(ntex)
    tmax = np.max(ntex)

    tmesh = a.fusionObjects([amesh, atex], method='FusionTexSurfMethod')
    tmesh.releaseAppRef()
    palette = self.palette
    atex.setPalette(palette, minVal=0., maxVal=1.)
    win = a.createWindow('Axial')
    bg_val = self.background_value
    if bg_val < tmin:
        bg_val = tmin
    elif bg_val > tmax:
        bg_val = tmax
    bg = atex.palette().refPalette().at(
        int(round(bg_val / (tmax - tmin)
                  * atex.palette().refPalette().getSizeX())))
    bg = [float(x) / 255. for x in bg]
    win.windowConfig(cursor_visibility=0, light={'background': bg})
    win.addObjects(tmesh)

    width = self.width
    height = self.height
    bbox = amesh.boundingbox()
    if self.keep_aspect_ratio:
        ratio = (bbox[1][1] - bbox[0][1]) / (bbox[1][0] - bbox[0][0])
        if width != 0:
            height = int(np.ceil(width * ratio))
        else:
            width = int(np.ceil(height / ratio))

    #if not self.output_rgb_image:
        ## build a depth map instead of a RGB image
        #mesh = a.toAimsObject(amesh)
        #vert = mesh.vertex()
        #for i, v in enumerate(vert):
            #v[2] = ntex[i]
        #amesh.setChanged()
        #amesh.notifyObservers()
        #im = win.snapshotImage(width, height, 8)

    im = win.snapshotImage(width, height)

    dtype = ntex.dtype
    scale = (tmax - tmin) / 255.

    if self.output_rgb_image:
        aim = qt_backend.qimage_to_np(im)
        vol = aims.Volume(aim.shape[1], aim.shape[0], dtype='RGBA')
        context.write(aim.shape, vol.shape)
        # slow copy
        for y in range(aim.shape[0]):
            for x in range(aim.shape[1]):
                vol.setValue(aims.AimsRGBA(*[int(v) for v in aim[y, x]]), x, y)
    else:
        aim = (qt_backend.qimage_to_np(im).astype(np.float32)[:, :, 0] * scale
              + tmin).astype(dtype)
        context.write('shape:', aim.shape)
        aim = np.array(aim) # copy array
        vol = aims.Volume(aim)
        vol.header()['palette'] = {'palette': self.palette}

    vsx = (bbox[1][0] - bbox[0][0]) / vol.getSizeX()
    vsy = (bbox[1][1] - bbox[0][1]) / vol.getSizeY()
    vol.header()['voxel_size'] = [vsx, vsy, 1., 1.]
    aims.write(vol, self.image_2d.fullPath())

    #return [tmesh, win]


def execution(self, context):

    return mainThreadActions().call(self.execution_mainthread, context)

