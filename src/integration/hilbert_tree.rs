use std::collections::HashMap;
use super::hilbert_model::HilbertModel;
use super::handle::Handle;
use num::BigUint;
use num::bigint::ToBigUint;
use num::traits::{ Zero, One };
use std::ops::Shl;

// Structures in this file refer to points using several coordinate systems:
//   - Domain Space. Floating point coordinates in the problem domain in D dimensions.
//   - Hilbert Axes. Unsigned numbers where values for all D dimensions range from zero to (2^bits) - 1.
//   - Hilbert Index. A one dimensional, big integer that results from applying the Hilbert Transform to the Hilbert Axes.
//   - Normalized Index. Ratio of the Hilbert Index to the length of the Hilbert curve in the range [0,1) as a floating point number.
//   - Handle. Identifies a segment of the Hilbert curve that resulted from a series of bisections.
//
// Functions must be supplied:
//   - to convert domain space points to and from Hilbert Indixes (via Hilbert Axes).
//   - to compute the distance between two domain space points.
//
// Since all segments will be derived via bisections of other segments, the Hilbert Fraction
// will always be a negative integer power of two: 1, 1/2, 1/4, 1/8, ...
// Because of this simplifying assumption, this fraction will be defined as a number of Cuts.
//   Fraction = (1/2)^cuts


/// A segment of the Hilbert curve.
/// 
/// Segments define the position of the start of the segment in domain coordinates, 
/// value of the integrand at that point,
/// and a handle to the point at the end of the segment. 
#[derive(Clone,Debug)]
struct Segment {
    /// Start of segment in Problem Domain coordinate space.
    pub start_point: Vec<f64>,

    /// Value of integrand evaluated at the start_point.
    pub value: f64,

    /// Reference to a node succeeding this Segment in Hilbert curve order that contains
    /// a Segment definition, hence not a LeftBranch or LeftLeaf.
    pub end_handle: Handle
}

impl Segment {
    /// Create an empty Segment.
    pub fn new_empty() -> Self {
        Segment {
            start_point: Vec::new(),
            value: f64::NAN,
            end_handle: Handle::None
        }
    }
}

/// Measures the value of the integrand at the start point of a segment, 
/// the distance from the start to end point in Problem space, and the slope.
pub struct SegmentMeasure {
    /// Value of the integrand
    pub value: f64,

    /// Distance between start and end points for this segment in Problem space.
    pub distance: f64,

    /// Change in value from start to end of segment divided by distance.
    pub slope: f64
}

/// A node in a HilbertTree (a binary tree).
/// 
/// New values only appear in the tree as RightBranch or RightLeaf nodes,
/// because a split reuses the Segment values from the parent in the LeftBranch or LeftLeaf.
/// That is why left nodes point to an ancestor for their values.
/// The ancestor a left node points to may be several generations up the tree. 
#[derive(Clone,Debug)]
enum TreeNode {
    /// A binary tree begins as a single RootLeaf with no children.
    RootLeaf(Box<Segment>),

    /// This is the root of the binary tree after its first split. 
    RootBranch(Box<Segment>),

    /// A branch node with children that is the left child of its parent.
    /// Only a handle is needed because its fields may be inferred from other nodes.
    /// The Handle points to the nearest ancestor that contains a Segment: a RightBranch or the RootBranch.
    LeftBranch(Handle),

    /// A branch node with children that is the right child of its parent.
    /// It requires more information about the segment's values than a LeftBranch.
    RightBranch(Box<Segment>),

    /// A leaf node with no children that is the left child of its parent.
    LeftLeaf(Handle),

    /// A leaf node with no children that is the right child of its parent.
    RightLeaf(Box<Segment>),

    /// A special marker for the last point in the Hilbert Curve.
    Terminus(Box<Segment>),

    /// Not a node
    Empty
}

impl TreeNode {
    pub fn is_empty(&self) -> bool {
        match self {
            Self::Empty => true,
            _ => false
        }
    } 

    /// Is the TreeNode a leaf in the tree or not?
    pub fn is_leaf(&self) -> bool {
        match self {
            Self::RootLeaf(_) => true,
            Self::LeftLeaf(_) => true,
            Self::RightLeaf(_) => true,
            _ => false
        }
    } 
}

/// A binary tree of segments of a Hilbert curve that can be recursively divided
/// and randomly searched according to a discrete distributiuon in log N time.
/// During division, each smaller segment is formed by bisecting its parent segment.
/// 
/// At any given time, the tree need not be balanced or full and 
/// does not satisfy the HEAP or SHAPE properties of a classical heap.
/// 
/// Nodes near the top of the tree (levels 0 to 15) are stored in a dense array (the heap)
/// while nodes farther from the root node are stored in sparse storage in a HashMap. 
/// 
/// Leaf nodes may be searched for using a normalized Hilbert index, a float in the range [0,1).
/// Only one leaf can contain that index, and it is found using a binary search down the binary tree.
/// By this means, if the normalized Hilbert index is generated randomly on the interval [0,1), 
/// the selection of the segment is proportional to a discrete probability distribution.
/// As nodes are split, this discrete probability distribution is effectively updated. 
/// 
/// The standard algorithm for sampling from a discrete distribution is the Vose Alias method.
/// It is faster than this algorithm when dealing with drawing from an unchanging distribution
/// (constant versus log N time), but rebuilding the Alias structure after an update to the 
/// probabilities is costly, making it unsuitable in this application.
pub struct HilbertTree<M> where M: HilbertModel {
    /// Transforms points between coordinate systems and evaluates the integrand for a point.
    model: M,

    /// Dense storage of nodes with smaller Handles
    heap: Vec<TreeNode>,

    /// Sparse storage of nodes with larger handles
    map: HashMap<u32,TreeNode>,

    /// Node for the last point along the curve.
    terminus: TreeNode

}

impl<M> HilbertTree<M> where M: HilbertModel {

    // ///////////////////////////////////
    //                                  //
    //   Construction & Initialization  //
    //                                  //
    // ///////////////////////////////////

    /// Construct a HilbertTree
    pub fn new(model: M) -> Self {
        const INIT:TreeNode = TreeNode::Empty;
        const HEAP_SIZE: usize = 65536;
        let mut tree = HilbertTree {
            model: model,
            heap: Vec::with_capacity(HEAP_SIZE),
            map: HashMap::new(),
            terminus: TreeNode::Terminus(Box::new(
                Segment::new_empty()
            ))
        };
        for _i in 0..HEAP_SIZE {
            tree.heap.push(INIT);
        }
        tree.init_root();
        tree
    }

    /// Initialize the root and terminus.
    fn init_root(&mut self) {
        let origin_index = BigUint::zero();
        let origin_point = self.model.to_point(&origin_index);
        let value = self.model.evaluate(&origin_point);
        let tp = self.terminus_point();
        let tp_value = self.model.evaluate(&tp);
        self.terminus = TreeNode::Terminus(Box::new(
            Segment {
                start_point: tp,
                value: tp_value,
                end_handle: Handle::None
            })
        );

        self.heap[0] = TreeNode::RootLeaf(Box::new(
            Segment {
                start_point: origin_point,
                value: value,
                end_handle: Handle::Terminus
            })
        );
    }

    fn terminus_point(&self) -> Vec<f64> {
        let shift = self.model.bits() as usize * self.model.dimensions() as usize;
        let terminus = BigUint::one().shl(shift) - BigUint::one();
        let terminus_point = self.model.to_point(&terminus);
        terminus_point
    }



    // //////////////////////////////////////////////
    //                                             //
    //                  Querying                   //
    //                                             //
    //  is_empty                                   //
    //  is_leaf                                    //
    //  get_value                                  //
    //  get_weight                                 //
    //  get_weighted_value                         //
    //  measure                                    //
    //  get_bisection_change                       //
    //  find_leaf                                  //
    //                                             //
    // //////////////////////////////////////////////

    /// Is the given node empty or does the handle not refer to a node at all?
    pub fn is_empty(&self, handle: Handle) -> bool {
        match handle {
            Handle::None => true,
            Handle::Terminus => false,
            Handle::Heap(i) => self.heap[i as usize].is_empty(),
            Handle::Map(i) => {
                match self.map.get(&i) {
                    None => true,
                    Some(node) => node.is_empty()
                }
            }
        }
    }

    /// Is the given node a leaf node?
    pub fn is_leaf(&self, handle: Handle) -> bool {
        match handle {
            Handle::None => false,
            Handle::Terminus => false,
            Handle::Heap(i) => self.heap[i as usize].is_leaf(),
            Handle::Map(i) => {
                match self.map.get(&i) {
                    None => false,
                    Some(node) => node.is_leaf()
                }
            }
        }
    }

    /// Get the value at the start point for the node.
    pub fn get_value(&self, handle: Handle) -> f64 {
        match self.get_value_node(handle) {
            TreeNode::RightBranch(segment) => segment.value,
            TreeNode::RightLeaf(segment) => segment.value,
            TreeNode::RootLeaf(segment) => segment.value,
            TreeNode::RootBranch(segment) => segment.value,
            TreeNode::Terminus(segment) => segment.value,
            _ => f64::NAN
        }
    }

    /// Get the weight for the segment in the inclusive range [0..1],
    /// or NAN if the Handle is empty. 
    pub fn get_weight(&self, handle: Handle) -> f64 {
        match handle.level() {
            Some(level) => 2.0_f64.powi(level as i32),
            None => f64::NAN
        }
    }

    /// Get the weighted value for the indicated segment, 
    /// which will be NAN for the terminus, an empty node, or a start point 
    /// for which the integrand function may not be evaluated.
    pub fn get_weighted_value(&self, handle: Handle) -> f64 {
        self.get_value(handle) * self.get_weight(handle)
    }

    /// Measure the value, length and slope of the segment.
    pub fn measure(&self, handle: Handle) -> SegmentMeasure {
        let start_handle = self.get_value_handle(handle);
        let start_segment = self.get_segment(start_handle);
        let start_point = &start_segment.start_point;
        let start_value = start_segment.value;

        let end_handle = self.get_end_node_handle(handle);
        let end_segment = self.get_segment(end_handle);
        let end_point = &end_segment.start_point;
        let end_value = end_segment.value;

        let distance = self.model.distance(start_point, end_point);
        SegmentMeasure {
            value: start_value,
            distance: distance,
            slope: (end_value - start_value) / distance
        }
    }

    /// Compute the change to the weighted sum of values attributable to the bisection
    /// of the indicated node into two shorter segments.
    /// If the handle refers to an empty or an as yet unsplit node, return None. 
    /// Weighted values that are NAN or infinity will be treated as zero.
    pub fn get_bisection_change(&self, handle: Handle) -> Option<f64> {
        if self.is_empty(handle) { return None; }
        let left_child = handle.get_left_child().unwrap();
        // Assume that if the left child is not empty, the right isn't either.
        if self.is_empty(left_child) { return None; } 
        let right_child = handle.get_right_child().unwrap();

        let parent_weighted_value = Self::nz(self.get_weighted_value(handle));
        let left_weighted_value = Self::nz(self.get_weighted_value(left_child));
        let right_weighted_value = Self::nz(self.get_weighted_value(right_child));

        Some(left_weighted_value + right_weighted_value - parent_weighted_value)
    }

    /// Walk the tree to find the handle to the leaf node whose range of values 
    /// includes the given normalized index. 
    /// 
    /// A normalized index must be in the range [0,1). 
    /// It indicates what fraction of the distance from the origin to 
    /// the final point on the Hilbert curve to travel.
    /// 
    /// By choosing a random or pseudo-random value for the normalized_index,
    /// you can sample this tree according to a discrete distribution 
    /// where each segment's probability of being chosen is proportional
    /// to its Hilbert curve length.
    pub fn find_leaf(&self, normalized_index: f64) -> Handle {
        if normalized_index < 0.0 || normalized_index >= 1.0 { return Handle::None; }

        // The root includes the normalized_index by default (if in range).
        let mut current_handle = Handle::new(0);  
        while !self.is_leaf(current_handle) {
            let left_child = current_handle.get_left_child().unwrap();
            if left_child.includes(normalized_index) {
                current_handle = left_child;
            }
            else {
                current_handle = current_handle.get_right_child().unwrap();
            }
        }
        current_handle
    }

    // //////////////////////////////////////////////
    //                                             //
    //            Querying Internals               //
    //                                             //
    //  get_node, get_value_node                   //
    //  get_value_handle, get_end_node_handle      //
    //  get_segment                                //
    //                                             //
    // //////////////////////////////////////////////

    /// Get a clone of the TreeNode referred to by the given Handle.
    fn get_node(&self, handle: Handle) -> TreeNode {
        match handle {
            Handle::None => TreeNode::Empty,
            Handle::Terminus => self.terminus.clone(),
            Handle::Heap(i) => {
                self.heap[i as usize].clone()
            },
            Handle::Map(i) => {
                match self.map.get(&i) {
                    None => TreeNode::Empty,
                    Some(node) => node.clone()
                }
            }
        }
    }

    /// Get a clone of the TreeNode referred to by the given Handle,
    /// or one of its ancestors, if the node is a Left node that
    /// only indirectly knows its value. 
    fn get_value_node(&self, handle: Handle) -> TreeNode {
        let node = self.get_node(handle);
        match node {
            TreeNode::LeftBranch(ancestor) => self.get_node(ancestor),
            TreeNode::LeftLeaf(ancestor) => self.get_node(ancestor),
            _ => node
        }
    }

    /// Get the Handle to the node that holds the point and value for the node indicated by this Handle,
    /// which may be the same node or an ancestor.
    fn get_value_handle(&self, handle: Handle) -> Handle {
        let empty = TreeNode::Empty;
        let node = match handle {
            Handle::None => &empty,
            Handle::Terminus => &self.terminus,
            Handle::Heap(i) => &self.heap[i as usize],
            Handle::Map(i) => {
                match self.map.get(&i) {
                    None => &empty,
                    Some(node) => node
                }
            }
        };
        match node {
            TreeNode::LeftBranch(ancestor) => *ancestor,
            TreeNode::LeftLeaf(ancestor) => *ancestor,
            _ => handle
        }
    }

    /// Find the node at the end of the segment indicated by the given Handle.
    fn get_end_node_handle(&self, handle: Handle) -> Handle {
        let empty = TreeNode::Empty;
        let handle_index = handle.get_index();
        let start_node = match handle {
            Handle::None => &empty,
            Handle::Terminus => &empty,
            Handle::Heap(i) => &self.heap[i as usize],
            Handle::Map(i) => {
                match self.map.get(&i) {
                    None => &empty,
                    Some(node) => &node
                }
            }
        };
        match start_node {
            TreeNode::LeftBranch(_) => Handle::new(handle_index.unwrap() + 1),
            TreeNode::LeftLeaf(_) => Handle::new(handle_index.unwrap() + 1),
            TreeNode::RightBranch(segment) => segment.end_handle,
            TreeNode::RightLeaf(segment) => segment.end_handle,
            TreeNode::RootBranch(_) => Handle::Terminus,
            TreeNode::RootLeaf(_) => Handle::Terminus,
            _ => Handle::None
        }
    }

    fn get_segment(&self, handle: Handle) -> Segment {
        match self.get_value_node(handle) {
            TreeNode::RightBranch(segment) => (*segment).clone(),
            TreeNode::RootLeaf(segment) => (*segment).clone(),
            TreeNode::RootBranch(segment) => (*segment).clone(),
            TreeNode::Terminus(segment) => (*segment).clone(),
            _ => Segment::new_empty()
        }
    }

    // ///////////////////////////////////
    //                                  //
    //           Modifying              //
    //                                  //
    //  bisect                          //
    //  bisect_index                    //
    //  set_node                        //
    //                                  //
    // ///////////////////////////////////

    /// Bisect the indicated node if it is a leaf and return a tuple with the Handles of the 
    /// left and right children and the change to the weighted sum of values that this split
    /// causes.
    pub fn bisect(&mut self, bisect_at: Handle) -> Option<(Handle, Handle, f64)> {
        /* 
            Steps to take: 
                1. Replace leaf being bisected with the appropriate branch
                     RootLeaf -> RootBranch, 
                     LeftLeaf -> LeftBranch, 
                     RightLeaf -> RightBranch
                2. Create right child, find the Hilbert point, convert to Problem space, evaluate integrand
                3. Set right child into the tree
                4. Create left child, reuse point and integrand of parent
                5. Set left child into the tree
                6. Compute change in integrand sum
        
        */
        if !self.is_leaf(bisect_at) {
            return None;
        } 

        let left_handle = bisect_at.get_left_child().unwrap();
        let right_handle = bisect_at.get_right_child().unwrap();
        let (point, value) = self.eval(right_handle).unwrap();

        match self.get_node(bisect_at) {
            TreeNode::LeftLeaf(value_handle) => {
                self.set_node(bisect_at, TreeNode::LeftBranch(value_handle));
                self.set_node(left_handle, TreeNode::LeftLeaf(value_handle));
                self.set_node(right_handle, TreeNode::RightLeaf(Box::new(
                    Segment {
                        start_point: point,
                        value: value,
                        end_handle: self.get_end_node_handle(bisect_at)
                    }
                )));
            },
            TreeNode::RightLeaf(segment) => {
                self.set_node(bisect_at, TreeNode::RightBranch(segment));
                self.set_node(left_handle, TreeNode::LeftLeaf(self.get_value_handle(bisect_at)));
                self.set_node(right_handle, TreeNode::RightLeaf(Box::new(
                    Segment {
                        start_point: point,
                        value: value,
                        end_handle: self.get_end_node_handle(bisect_at)
                    }
                )));
            },
            TreeNode::RootLeaf(segment) => {
                let root_end_handle = segment.end_handle;
                self.set_node(bisect_at, TreeNode::RootBranch(segment));
                self.set_node(left_handle, TreeNode::LeftLeaf(bisect_at));
                self.set_node(right_handle, TreeNode::RightLeaf(Box::new(
                    Segment {
                        start_point: point,
                        value: value,
                        end_handle: root_end_handle
                    }
                )));
            },
            _ => { return None; }
        };
        Some((left_handle, right_handle, self.get_bisection_change(bisect_at).unwrap()))
    }

    pub fn bisect_index(&mut self, normalized_index: f64) -> Option<(Handle, Handle, f64)>  {
        let node_to_bisect = self.find_leaf(normalized_index);
        self.bisect(node_to_bisect)
    }

    fn set_node(&mut self, handle: Handle, node: TreeNode) -> bool {
        match handle {
            Handle::None => false,
            Handle::Terminus => { self.terminus = node; true },
            Handle::Heap(i) => { 
                self.heap[i as usize] = node; true 
            },
            Handle::Map(i) => { self.map.insert(i, node); true } 
        }
    }

    // ///////////////////////////////////
    //                                  //
    //           Utility                //
    //                                  //
    //  nz                              //
    //  eval                            //
    //  get_hilbert_index               //
    //  draw                            //
    //                                  //
    // ///////////////////////////////////

    /// Convert NAN or infinity into zero.
    fn nz(x: f64) -> f64 {
        if x.is_finite() { x }
        else { 0.0 }
    }

    /// Find the point in the Problem coordinate space for the given Handle
    /// and evaluate the integrand function to get the value. 
    fn eval(&self, node: Handle) -> Option<(Vec<f64>, f64)> {
        match self.get_hilbert_index(node) {
            Some(hilbert_index) => {
                let problem_space = self.model.to_point(&hilbert_index);
                let value = self.model.evaluate(&problem_space);
                Some((problem_space, value))
            },
            None => None
        }
    }

    /// Translate the handle into a HilbertIndex. 
    /// 
    /// ```text
    ///                        bits x dimensions
    ///    Curve Length  L = 2
    ///    
    ///                                        -level
    ///    Hilbert Index I = L x numerator x 2
    ///         
    ///                                    bits x dimensions - level
    ///                    = numerator x 2
    /// ```
    fn get_hilbert_index(&self, handle: Handle) -> Option<BigUint> {
        match (handle.numerator(), handle.level()) {
            (Some(numerator), Some(level)) => {
                let power = self.model.bits() as usize * self.model.dimensions() as usize - level as usize;
                let big_numerator: BigUint = numerator.to_biguint().unwrap();
                Some(big_numerator.shl(power))
            },
            _ => None
        }
    }


    fn node_name(node: Handle) -> String {
        let level = node.level().unwrap() as usize;
        let numerator = node.numerator().unwrap();
        let denominator = 1_usize << level;   
        let name = if level == 0 {
            format!("[0,1)")
        }
        else {
            format!("[{}/{},{}/{})", numerator, denominator, numerator + 1, denominator)
        };
        name
    }

    /// Format the tree as a DOT string that Graphviz can render as a Binary tree diagram.
    pub fn draw(&self, node: Handle, max_depth: usize, decimals: usize, output: &mut String) {
        if self.is_empty(node) { return; }
        let current_level = node.level().unwrap() as usize;
        if current_level > max_depth { return; }

        let value = self.get_value(node);
        let name = Self::node_name(node);
        let left_child = node.get_left_child().unwrap();
        let right_child = node.get_right_child().unwrap();
        let label = format!("{:.2$}\\n{}", value, name, decimals);   

        if current_level == 0 {
            // The header for the DOT.
            output.push_str(&format!("
digraph BST {{
    node [fontname=\"Arial\"];

    \"{}\" [shape=doublecircle,label=\"{}\"]", name, label));
        }
        else { 
            if self.is_leaf(node) {
                output.push_str(&format!("    \"{}\" [shape=octagon,label=\"{}\"];\n", name, label));
            }
            else {
                output.push_str(&format!("\n    \"{}\" [label=\"{}\"];", name, label));
            }
        }
        
        if !self.is_leaf(node) {
            // Links to child nodes.
            output.push_str(&format!("
    \"{}\" -> \"{}\";
    \"{}\" -> \"{}\";
", name, Self::node_name(left_child), name, Self::node_name(right_child)));

            // Recurse to child nodes.
            self.draw(left_child, max_depth, decimals, output);
            self.draw(right_child, max_depth, decimals, output);
        }

        // The footer for the DOT.
        if current_level == 0 {
            output.push_str("}\n");
        }
    }

}

#[cfg(test)]
mod tests {
    use super::Handle;
    use super::HilbertTree;
    use super::super::hilbert_model::{HilbertModel,Quantizer,DimensionRange};

    #[test]
    fn handle_level() {
        assert_eq!(Handle::new(0).level().unwrap(), 0);
        assert_eq!(Handle::new(1).level().unwrap(), 1);
        assert_eq!(Handle::new(2).level().unwrap(), 1);
        assert_eq!(Handle::new(3).level().unwrap(), 2);
        assert_eq!(Handle::new(4).level().unwrap(), 2);
        assert_eq!(Handle::new(5).level().unwrap(), 2);
        assert_eq!(Handle::new(6).level().unwrap(), 2);
        assert_eq!(Handle::new(u16::MAX as u32).level().unwrap(), 16);
    }

    #[test]
    fn handle_numerator() {
        assert_eq!(Handle::new(0).numerator().unwrap(), 0);
        assert_eq!(Handle::new(1).numerator().unwrap(), 0);
        assert_eq!(Handle::new(2).numerator().unwrap(), 1);
        assert_eq!(Handle::new(3).numerator().unwrap(), 0);
        assert_eq!(Handle::new(4).numerator().unwrap(), 1);
        assert_eq!(Handle::new(5).numerator().unwrap(), 2);
        assert_eq!(Handle::new(6).numerator().unwrap(), 3);
    }


    #[test]
    fn handle_includes() {
        assert!(Handle::new(0).includes(0.25));
        assert!(Handle::new(1).includes(0.25));
        assert!(Handle::new(2).includes(0.75));
        assert!(Handle::new(3).includes(0.125));
        assert!(Handle::new(4).includes(0.25));
        assert_eq!(Handle::new(5).includes(0.25), false);
        assert_eq!(Handle::new(6).includes(f64::NAN), false);
    }

    #[test]
    fn hilbert_tree_bisect_index() {
        // Curious about this text?
        // Go to https://dreampuf.github.io/GraphvizOnline and then next
        // paste it in and you will see 
        // a nicely labeled binary tree.
        let expected_graph = r###"
digraph BST {
    node [fontname="Arial"];

    "[0,1)" [shape=doublecircle,label="1.000\n[0,1)"]
    "[0,1)" -> "[0/2,1/2)";
    "[0,1)" -> "[1/2,2/2)";

    "[0/2,1/2)" [label="1.000\n[0/2,1/2)"];
    "[0/2,1/2)" -> "[0/4,1/4)";
    "[0/2,1/2)" -> "[1/4,2/4)";
    "[0/4,1/4)" [shape=octagon,label="1.000\n[0/4,1/4)"];
    "[1/4,2/4)" [shape=octagon,label="1.000\n[1/4,2/4)"];

    "[1/2,2/2)" [label="262145.000\n[1/2,2/2)"];
    "[1/2,2/2)" -> "[2/4,3/4)";
    "[1/2,2/2)" -> "[3/4,4/4)";
    "[2/4,3/4)" [shape=octagon,label="262145.000\n[2/4,3/4)"];

    "[3/4,4/4)" [label="524287.500\n[3/4,4/4)"];
    "[3/4,4/4)" -> "[6/8,7/8)";
    "[3/4,4/4)" -> "[7/8,8/8)";
    "[6/8,7/8)" [shape=octagon,label="524287.500\n[6/8,7/8)"];

    "[7/8,8/8)" [label="196608.000\n[7/8,8/8)"];
    "[7/8,8/8)" -> "[14/16,15/16)";
    "[7/8,8/8)" -> "[15/16,16/16)";
    "[14/16,15/16)" [shape=octagon,label="196608.000\n[14/16,15/16)"];
    "[15/16,16/16)" [shape=octagon,label="1.000\n[15/16,16/16)"];
}
"###;
        let model = create_test_model();
        let mut tree = HilbertTree::new(*model);

        tree.bisect_index(0.25);
        tree.bisect_index(0.75);
        tree.bisect_index(0.125);
        tree.bisect_index(0.875);
        tree.bisect_index(0.9);
        let mut graph = String::new();
        tree.draw(Handle::new(0), 6, 3, &mut graph);
        println!("Graphviz DOT:\n{}", graph);
        assert_eq!(expected_graph, &graph);
    }


    fn create_test_model() -> Box<impl HilbertModel> {
        let f = |p: &[f64]| -> f64 { 
            let y = p[0]*p[1] + 1.0;
            // println!("eval f(p) for [{},{}] = {}", p[0], p[1], y);
            y
        };
        let ranges = vec![
            DimensionRange::new(0.0, 1024.0, 1.0, 20),
            DimensionRange::new(0.0, 1024.0, 1.0, 20),
        ];
        Box::new(Quantizer::new(ranges, 20, f))
    }
}

