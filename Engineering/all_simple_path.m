function p = all_simple_path(Edges, Ntot, Source, Sink, CUTOFF)
% Find all simple paths between the source and sink node. The list of all
% simple paths is returned as an cell array "p";
%
% Input:
%
%   Edges = [Node_s Node_e]; N x 2 vector of the starting and end nodes
% 
%   Ntot; total number of nodes in the graph (start from node 1 by default)
%         if Ntot is empty, then assign the maximum node index to Ntot
%
%   Source, Sink; Index of the source and sink node in the graph
%
%   CUTOFF; only find path whose length is shorter than CUTOFF
%           default CUTOFF = min(10,Ntot)

Node_s = Edges(:,1);
Node_e = Edges(:,2);
if isempty(Ntot); Ntot = max(Node_e(:)); end
if ~exist('CUTOFF','var') || isempty(CUTOFF); CUTOFF = min(10,Ntot); end



global Stack Layer Visited Paths PathCount
Paths = {};
Stack   = [];
Visited = zeros(Ntot,1);
Stack(1) = Source;
Visited(Source) = 1;
Layer     = 1;
PathCount = 0;

while Layer > 0 
    
    DFSutil(Node_s, Node_e, Source, Sink, CUTOFF);
    
end

p = Paths;


function DFSutil(Node_s, Node_e, v, Sink, CUTOFF)
global Stack Layer Visited Paths PathCount
if v == Sink 
%     disp(Stack);
    PathCount = PathCount + 1;
    Paths{PathCount} = Stack;
    Visited(Stack(Layer)) = 0;
    Stack(Layer) = [];
    Layer = Layer - 1;
    return
end

if Layer >= CUTOFF
    Visited(Stack(Layer)) = 0;
    Stack(Layer) = [];
    Layer = Layer - 1;
    return
end

% find all connected nodes
v1 = Node_s(Node_e == v);
v2 = Node_e(Node_s == v);
vn = [v1; v2];
for i = 1:numel(vn)
    if Visited(vn(i)) == 0
        Layer = Layer + 1;
        Stack(Layer) = vn(i);
        Visited(vn(i))  = 1;
%         disp(Stack)
        DFSutil(Node_s, Node_e, vn(i), Sink , CUTOFF)
    end
end
Visited(Stack(Layer)) = 0;
Stack(Layer) = [];
% disp(Stack);
Layer = Layer - 1;
    
function printStack()
global Stack
for i = 1:numel(Stack)
    fprintf('%d ',Stack(i))
end
fprintf('\n')
