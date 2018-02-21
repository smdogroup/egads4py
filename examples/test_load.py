from egads4py import egads

ctx = egads.context()
ctx.setOutLevel(2)
ego = ctx.loadModel('misc1.step')

geo, oclass, mtype, lim, body, senses = ego.getTopology()

for b in body:
    faces = b.getBodyTopos(egads.FACE)

    for f in faces:
        loops = b.getBodyTopos(egads.LOOP, ref=f)

        for l in loops:
            edges = b.getBodyTopos(egads.EDGE, ref=l)

            print edges



    
