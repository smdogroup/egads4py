
class egadsInterNode
{
public:
  TopoDS_Vertex        newNode;         // the new Edge Node
  TopoDS_Vertex        Node[2];         // hit existing Node src/tool
  TopoDS_Edge          Edge[2];         // position is at t on this Edge
  double               t[2];            // t for splitting Edge src/tool
};


class egadsInterEdge
{
public:
  TopoDS_Edge          newEdge;         // the new Edge
  TopoDS_Face          oFace[2];        // the original Face cut by EDGE
  TopoDS_Face          mFace[2];        // the Face cut by EDGE in src/tool
  TopoDS_Edge          Edge[2];         // coincident Edge from src/tool
  egadsInterNode       start;           // the beginning Edge condition
  egadsInterNode       end;             // the ending Edge condition
  int                  prev;            // previous segment index (1 bias)
  int                  next;            // next segment index (0 is open)
};



extern int EG_intersect(TopoDS_Shape src, TopoDS_Shape tool, int outLevel,
                        int *nInter, egadsInterEdge **inters);

