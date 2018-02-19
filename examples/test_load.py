from egads4py import egads

ego = egads.pyego()
ego.loadModel('misc1.step')

geo, oclass, mtype, lim, children, senses = ego.getTopology()

for c in children:
    f = c.getBodyTopos(egads.FACE)
    print len(f)

print ego

    
