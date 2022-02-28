/// Handle to a node in the HilbertTree, a binary tree.
/// 
/// Because the classic heap layout is used, pointers are not necessary to connect
/// parents to children and vice versa. The Handle to a parent may be directly computed 
/// from its child's integer index and vice versa using simple integer arithmetic.
/// 
/// Because we are using a tree to represent bisections of a line segment into children of equal lengths,
/// instead of arbitrary divisions, the Handle may also be used to derive the Normalized Hilbert Index,
/// using arithmetic based solely on the Handle's integer index. By making the pointers between 
/// tree nodes implicit, not explicit, memory use is reduced and tricky Borrows and lifetimes are not needed.
#[derive(Clone,Copy,Debug)]
pub enum Handle {
    /// A Heap handle is an index to a value stored in the dense array, for indices in [0,65535).
    /// Array will not have a Slot 65535, as it begins the next row of the heap.
    Heap(u16),

    /// A Map handle is the key to a value stored in the sparse Map, for indices >= 65536. 
    Map(u32),

    /// A Terminus handle refers to the maximum value of the Hilbert Index, the one at the end of the Hilbert curve.
    Terminus,

    /// Handle pointing to nothing (used only for the segment Terminus)
    None
} 

impl Handle {
    pub fn new(index: u32) -> Handle {
        if (index + 1) < 1_u32 << 16 {
            Handle::Heap(index as u16)
        }
        else {
            Handle::Map(index)
        }
    }

    pub fn get_index(&self) -> Option<u32> {
        match self {
            Handle::Heap(i) => Some(*i as u32),
            Handle::Map(i) => Some(*i),
            _ => None
        }
    }

    pub fn increment(&self) -> Handle {
        match self.get_index() {
            Some(index) => Handle::new(index + 1),
            None => Handle::None
        }
    }

    /// Get the level in the binary tree for this Handle. 
    /// Level zero is the root node, level one the root's two children, level two is level one's four children, etc.
    pub fn level(&self) -> Option<u16> {
        // Since the tree is arranged as a heap, the floor(log2(index+1)) gives us the level. 
        // Use leading_zeros() function because it uses a low-level intrinsic function that is faster than system log combined with floor.
        match self {
            Handle::Heap(0) => Some(0),
            Handle::Heap(i) => Some((15 - (i+1).leading_zeros()) as u16),
            Handle::Map(i) => Some((31 - (i+1).leading_zeros()) as u16),
            _ => None
        }
    }

    pub fn get_left_child(&self) -> Option<Handle> {
        match self {
            Handle::Heap(i) => Some(Handle::new(1 + *i as u32 * 2)),
            Handle::Map(i) => Some(Handle::new(1 + *i * 2)),
            _ => None
        }
    }

    pub fn get_right_child(&self) -> Option<Handle> {
        match self {
            Handle::Heap(i) => Some(Handle::new(2 + *i as u32 * 2)),
            Handle::Map(i) => Some(Handle::new(2 + *i * 2)),
            _ => None
        }
    }

    /// Get the numerator of the fraction that locates the starting point of a 
    /// Hilbert curve segment proportional to its length. 
    /// The denominator is 2^level. 
    /// 
    /// To get the Hilbert Index I, if the length of the Hilbert curve is L: 
    /// 
    /// ```text
    ///                          -level
    ///    I = L * numerator * 2
    ///      
    /// ```
    /// Since L is always a power no smaller that level, I is always an integer.
    pub fn numerator(&self) -> Option<u32> {
        let level_opt = self.level();
        match (self, level_opt) {
            (_, None) => None,
            (Handle::Heap(0), Some(_level)) => Some(0),
            (Handle::Heap(i), Some(level)) => Some(*i as u32 + 1 - (1_u32 << level)),
            (Handle::Map(i), Some(level)) => Some(*i as u32 + 1 - (1_u32 << level)),
            _ => None
        }
    }  
    
    /// Return the inclusive low and exclusive high values as normalized Hilbert indices.
    /// 
    /// The value 0.0 corresponds to the start of the curve for Hilbert Axes [0,0,..,0]. 
    /// The value 1.0 corresponds to a point just beyond the end of the curve (hence is exclusive).
    /// (NAN,NAN) will be returned for None Handles.
    pub fn range(&self) -> (f64,f64) {
        match self.numerator() {
            None => (f64::NAN,f64::NAN),
            Some(numerator) => {
                let level = self.level().unwrap();
                let denominator = (1_u64 << level) as f64;
                (numerator as f64 / denominator, (numerator + 1) as f64 / denominator)
            }
        }
    }

    /// Decide if the normalized_index falls in the range defined by this Handle.
    pub fn includes(&self, normalized_index: f64) -> bool {
        if normalized_index < 0.0 || normalized_index >= 1.0 || !normalized_index.is_finite() {
            false
        }
        else {
            let (low_inclusive, high_exclusive) = self.range();
            normalized_index >= low_inclusive && normalized_index < high_exclusive            
        }
    }

}
