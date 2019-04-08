from egads4py import egads

# Create the egads context
ctx = egads.context()

# Load in a file, sew all the faces together then write the model
# back to a new file
filename = 'ucrm_9_surfaces.step'
crm = ctx.loadModel(filename)

all_faces = []
for body in crm.getChildren():
    all_faces.extend(body.getBodyTopos(egads.FACE))

model = ctx.sewFaces(all_faces, toler=1e-3, manifold=False)
model.saveModel('ucrm_9_model.step', overwrite=True)
