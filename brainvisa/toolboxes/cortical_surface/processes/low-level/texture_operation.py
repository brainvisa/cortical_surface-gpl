from brainvisa.processes import *

name = 'Texture operation'
userLevel = 0

signature = Signature(
    'textures', ListOf(ReadDiskItem('Texture', 'aims texture formats')),
    'operation', Choice('Formula', 'Min', 'Max'),
    'formula', String(),
    'output_texture', WriteDiskItem('Texture', 'aims texture formats'),
    'output_is_first_input', Boolean(),
)


def initialization(self):
    pass


def execution(self, context):
    from soma import aims
    import numpy as np

    T = []
    NT = []
    tex = {'T': T, 'NT': NT}
    for i, tname in enumerate(self.textures):
        T.append(aims.read(tname.fullPath()))
        tex['T%d' % i] = T[-1]
        try:
            ntex = np.asarray(T[-1][0])
        except:
            ntex = None
        NT.append(ntex)
        tex['NT%d' % i] = NT[-1]

    context.write(tex)

    if self.operation == 'Formula':
        if self.output_is_first_input:
            exec(self.formula, globals(), tex)
            res = T[0]
        else:
            res = eval(self.formula, globals(), tex)
    elif self.operation == 'Min':
        res = np.min(NT)
    elif self.operation == 'Max':
        res = np.max(NT)

    if type(res) is np.ndarray:
        otex = aims.TimeTexture(dtype=res.dtype)
        otex[0].assign(res)
        res = otex

    aims.write(res, self.output_texture.fullPath())

