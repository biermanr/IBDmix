March 17th 2026
---------------

I'm walking through IBDmix source code carefully now that I have
a slightly better understanding of (a) C++ (b) IBDmix itself and
(c) a desire to see if IBDmix can be parallelized (d) a desire
to see if IBDmix can use bcftools for random access.

I think I've done this at least once before, but as I look through it
now it all seems new, so I guess it didn't stick. I'll write my notes here.


So main.cc is what's used to create the `ibdmix` executable which is what
I'm most interested in. It starts by creating a `GenotypeReader` and then
creating an `IBD_collection`.

`IBD_collection::initialize` takes as input the `GenotypeReader` and creates
a vector of `IBD_segment`s called `IBDs`, one for each sample but not the archaic.
It does this by using the SampleMapper class which handles user-specified sample subsets.

Initializing an `IBD_segment` requires a "name", "threshold", "pool", and
"exclusive" flag. These are all simple except for "pool" which I'm still
trying to figure out.
- "name" is the sample name
- "threshold" is the LOD-threshold specified by args
- "exclusive" is a bool for using `(start, end)` or `(start, end]`
- "pool"  I'm still figuring out but it gets passed by reference
  to all the `IBD_segments` so I think it's a single shared obj

I'm surprised that `IBDs.emplace_back(sample, threshold, &pool, exclusive_end);`
in `IBD_Collection::initialize` is valid since it seems to be implicitly creating
`IBD_segment` objects and adding them to the IBDs vector in a single command. Didn't
know this was possible with `c++`.

The rest of main.cc proceeds by alternating calls to GenotypeReader::update()
and `IBD_collection::update()`.

`GenotypeReader::update()` reads the next line of the genotype file,
checks if the locus is within the mask file then calculates does more
processing with `GenotypeReader::process_line_buffer()`. This includes
finding the frequency at this locus using `Genotype_Reader::find_frequency()`
which uses only modern genotypes, not the archaic. Then it uses a `LODCalculator`,
which is a member of the `GenotypeReader` class, to calculate the lodscore of
this locus for each individual. Since there are only a handful of possible lodscores
per locus because given the allele frequency, it only depends on the modern/ancient genotypes.
That's why there is a `calculator.update_lod_cache` and then a `calculator.calculate_lod` because
IBDmix is caching the scores to avoid recalculating the same thing over and over.
These results are stored in the `lod_scores` member.

`IBD_collection::update()` takes as input the `GenotypeReader` and calls
`IBD_segment::add_lod()` on each of the `IBD_segments` in the `IBDs` vector.
This gets its info from the `GenotypeReader` `lod_scores` which were calculated
in the `GenotypeReader::update()` function.

This call to `IBD_segment::add_lod()` calls `IBD_sgement::add_node` which uses
the `IBD_pool->get_node` function to create(?) a new `IBD_node` and "push"
it onto the `IBD_stack` which is named `segment`. I think the `IBD_pool` is
being smart about memory and is re-using prior `IBD_Node` memory allocations, but I need
to look more into this. NOTE that if this is true, `IBD_pool` might make it complicated
to parallelize ibdmix, or `IBD_pool->get_node` would need to have atomic operations? It
looks like each `IBD_Segment` has it's own `IBD_pool`.

Any, this `IBD_Node` get's pushed onto `IBD_Segment`'s `segment`
member, which is an `IBD_Stack`.  The `IBD_stack` is a singly-linked list of
`IBD_Node` which is a data struct:
```
struct IBD_Node {
  double cumulative_lod, lod;
  uint64_t position;
  unsigned char bitmask;
  IBD_Node *next;
};
```
Note that `IBD_Node` does not need a `chrom` field because all nodes
in the same stack are guaranteed to be on the same `chrom`.

During the `IBD_segment::add_node`, a new node of the current position
becomes the new head of the singly-linked list of `IBD_Node`s that make
up the `IBD_stack`. The `cumulative_lod` of the new node is the sum of
it's own lod and the lod of the prior top.

After the `IBD_node` is pushed to `IBD_stack`, there are some logically
cases to consider such as:
* We only have one node on the stack, so we're starting a new segment
  -> We have to re-initialize the recorders to set all their values to 0
  -> Then record this new locus for all recorders

* `segment.topIsNewMax()` the best cumulative lod score is now the just-added node.
  -> NOTE I'm not very sure about these descriptions
  -> We set the end of the segment to be the new top which means that
     we've expanded the length of the segment because this is a "good" node
     that we just added.
  -> Then we run `pool->reclaim_segment(&segment)` but I'm not sure why.
     This calls `IBD_Stack::getSegmentFrom(IBD_Stack *other)` which I'm
     getting confused about. Seems like it's joining the pool stack and the segment stack?

* `segment.reachedMax()` means that by adding this last node, the cumulative lod
   at the top is now negative.
   -> We write-out the best segment from start to end
   -> Then do something with unprocessing the segment which involves
      reversing the linked list
   -> Then we `pool->reclaim_stack(&segment)` which is somehow returning nodes
      back to the pool?
   -> Finally it looks like we are cleaning up the unprocessed segment?

Side note that the `pool` member of `IBD_Pool` is actually an `IBD_Stack` itself,
which I guess makes sense since I think `IBD_Pool` is trying to keep track of
already allocated nodes for reuse?

Main take-aways:
> Main loop alternates between (a) fetching the next line of the genotype file,
  calculating the lod per sample and (b) adding these lods to IBD_collection.
> A single `IBD_collection` is created which has a vector of `IBD_segments`,
> one for each sample. These `IBD_segments` each have an `IBD_stack` which is
> a singly-linked list of `IBD_node`s.

Things I still don't understand
- The purpose of the `IBD_pool`
- How the `IBD_stack` is creating the introgressed regions


May 6th 2026
-------------

Filling in the gaps from above, specifically around IBD_Pool, IBD_Stack,
and how introgressed regions actually get detected and emitted.


### Correction: all segments share ONE pool

Earlier I wrote "It looks like each IBD_Segment has its own IBD_pool" but
that's wrong. `IBD_Collection` owns a single `IBD_Pool pool` member
(IBD_Collection.h:25), and passes `&pool` to every `IBD_Segment` it
creates. So all samples share one pool. This is the key performance
trick -- node memory is recycled globally.

For parallelization, this means the pool IS a bottleneck/shared resource.
Either each thread would need its own pool, or get_node/reclaim would
need synchronization.


### What IBD_Pool actually is

IBD_Pool is a slab allocator / free-list for IBD_Node objects. Its purpose
is to avoid calling malloc/free thousands of times per sample per chromosome.
Instead, it pre-allocates a contiguous array of nodes and hands them out
one at a time.

The pool's internal state is just:
- `IBD_Stack pool` -- a free-list of available nodes, stored as a linked list
- `vector<IBD_Node*> alloc_ptrs` -- tracks the raw allocations for cleanup
- `int buffer_size` -- how many nodes to allocate next (doubles each time)

On construction, it malloc's 1024 nodes, links them into a chain via their
`next` pointers, and pushes them onto the `pool` stack. That's the initial
free-list.

```
pool (IBD_Stack):
  top -> [node0] -> [node1] -> [node2] -> ... -> [node1023] -> nullptr
```

#### get_node(position, lod, bitmask)

Pops one node off the free-list, fills in its fields, returns it.
If the free-list is empty, it malloc's another batch (doubling in size
each time: 1024, 2048, 4096, ...) and links those into the free-list
before popping.

```
BEFORE get_node:
  pool: top -> [A] -> [B] -> [C] -> ...

AFTER get_node returns A:
  pool: top -> [B] -> [C] -> ...
  returned: [A] with position/lod/bitmask filled in, next=nullptr
```

#### reclaim_node(node)

Pushes a single node back onto the free-list. Used when we immediately
reject a node (e.g., negative LOD on an empty segment).

#### reclaim_stack(stack)

Moves ALL nodes from the given stack to the free-list. Used when a
segment is fully emitted and all its nodes are done.

```
BEFORE reclaim_stack:
  pool:    top -> [X] -> [Y] -> ...
  segment: top -> [A] -> [B] -> [C]
                                 ^start
AFTER (getAllFrom):
  pool:    top -> [A] -> [B] -> [C] -> [X] -> [Y] -> ...
  segment: empty (top=start=end=nullptr)
```

The trick is that getAllFrom splices the segment's linked list directly
onto the pool's list -- it just repoints a few pointers, O(1) work.
It attaches the segment's start->next to pool's old top, then sets
pool's top to segment's top.

#### reclaim_segment(stack)

This is the confusing one. It reclaims only the INTERIOR nodes of the
stack -- nodes between `end` and `start`, exclusive of both. These are
"used up" nodes that we've already accumulated past and don't need.

```
BEFORE reclaim_segment:
  segment: top/end -> [300] -> [200] -> [100]
                                         ^start
  (nodes to reclaim: [200], between end and start)

  pool: top -> [X] -> [Y] -> ...

AFTER (getSegmentFrom):
  segment: top/end -> [300] -> [100]
                                ^start
  pool:    top -> [200] -> [X] -> [Y] -> ...
```

Why? Because once end advances to a new position, we know the segment
stretches from start to end. The intermediate nodes only existed to track
the running sum. Their LOD values have already been accumulated into end's
cumulative_lod. We keep start (need its position) and end (need its
position + cumulative LOD), but everything in between is garbage we can
recycle.


### How introgressed regions are detected (the max-subarray algorithm)

The segment detection is an adaptation of Kadane's maximum subarray
algorithm. For each sample, we're looking for contiguous stretches of
genomic positions where LOD scores accumulate to a high positive value --
these are the introgressed regions.

The core idea:
1. Push LOD scores onto the stack, accumulating a running sum
2. Track where the running sum hits its maximum
3. When the running sum goes NEGATIVE, the current "segment" is done
4. If the max cumulative LOD exceeded the threshold, emit it
5. Reprocess the "tail" nodes after the maximum (they might start a new segment)

Here's a complete worked example. Suppose we have sites at positions
100-500 with LOD scores +2, +1, +3, -1, -8, and threshold=3.

#### Step 1: Push pos=100, lod=+2

```
segment: top -> [100: lod=+2, cum=2]
                 ^start ^end
```

isSingleton -> initialize recorders, record this node.
This is the start of a potential segment.

#### Step 2: Push pos=200, lod=+1

```
segment: top -> [200: lod=+1, cum=3] -> [100: lod=+2, cum=2]
                                          ^start ^end
```

Note the stack is in REVERSE genomic order (newest on top).
cum_lod of 200 = prev_cum(2) + lod(1) = 3.

topIsNewMax? cum(3) >= end_cum(2) -> YES!
  - Record stats for node 200
  - setEnd: end moves to top (node 200)
  - reclaim_segment: reclaim nodes between end(200) and start(100)
    -> they're adjacent, nothing to reclaim

```
segment: top/end -> [200: cum=3] -> [100]
                                     ^start
```

#### Step 3: Push pos=300, lod=+3

```
segment: top -> [300: lod=+3, cum=6] -> [200: cum=3] -> [100]
                                          ^end            ^start
```

topIsNewMax? cum(6) >= end_cum(3) -> YES!
  - Record stats for node 300
  - setEnd: end = top (node 300)
  - reclaim_segment: reclaim nodes between end(300) and start(100)
    -> [200] is between them, gets returned to pool!

```
segment: top/end -> [300: cum=6] -> [100]     pool gained: [200]
                                     ^start
```

This is the key insight about reclaim_segment: node [200] did its job
(its LOD was folded into the cumulative sum), and since end has advanced
past it, we'll never need it again. Recycle it.

#### Step 4: Push pos=400, lod=-1

```
segment: top -> [400: lod=-1, cum=5] -> [300: cum=6] -> [100]
                                          ^end            ^start
```

topIsNewMax? cum(5) >= end_cum(6) -> NO (5 < 6)
reachedMax? cum(5) < 0 -> NO
Nothing happens. Node 400 just sits on top. It's dragging the
cumulative sum down but not enough to end the segment.

#### Step 5: Push pos=500, lod=-8

```
segment: top -> [500: lod=-8, cum=-3] -> [400: cum=5] -> [300: cum=6] -> [100]
                                                           ^end            ^start
```

topIsNewMax? cum(-3) >= end_cum(6) -> NO
reachedMax? cum(-3) < 0 -> YES! The running sum went negative.

This triggers the segment finalization:

**5a. Emit the segment:**
  endLod(6) >= threshold(3) -> YES, write output!
  Output: `sample_name  chrom  100  300  6.0`
  (start=100, end=300, LOD=6)

**5b. getUnprocessed() -- separate the tail from the emitted segment:**

First, reverse the stack so it's in genomic order:
```
BEFORE reverse:
  top -> [500] -> [400] -> [300] -> [100]
                            ^end     ^start

AFTER reverse:
  top -> [100] -> [300] -> [400] -> [500]
                   ^end              ^start
```

Then split at end: everything after end is "unprocessed"
```
  unprocessed: top -> [400] -> [500]
                                ^start

  segment:     top -> [100] -> [300]     (start=end=[300])
                                ^start/end
```

**5c. reclaim_stack(&segment):**
  Nodes [100] and [300] go back to pool. Segment is now empty.

**5d. Reprocess the unprocessed tail:**
  Pop [400] (lod=-1):
    segment is empty AND lod < 0 -> reclaim_node, skip.

  Pop [500] (lod=-8):
    segment is empty AND lod < 0 -> reclaim_node, skip.

All nodes are back in the pool, one segment was emitted. Done!


### What about adjacent segments?

The reprocessing step is what handles back-to-back introgressed regions.
If the tail after the maximum contained positive LODs, they'd start
building a NEW segment. For example, if pos=500 had lod=+4 instead of
-8, it would survive reprocessing and become the start of a fresh segment.


### Putting it all together: node lifecycle

```
                malloc'd slab
                     |
                     v
    IBD_Pool free-list: [available nodes linked together]
         |                          ^          ^
         | get_node()               |          |
         v                          |          |
    IBD_Segment.segment (IBD_Stack) |          |
         |                          |          |
         |  reclaim_segment()  -----+          |
         |  (interior nodes recycled           |
         |   when end advances)                |
         |                                     |
         |  reclaim_stack() -------------------+
         |  (all nodes recycled when
         |   segment is emitted or discarded)
         |
         |  reclaim_node() --------------------+
         |  (single node recycled when         |
         |   immediately rejected, e.g.        |
         |   negative lod on empty segment)    |
         |                                     |
         +-------------------------------------+
               nodes cycle back to pool
```

Nodes flow: pool -> segment stack -> pool -> segment stack -> ...
They're never freed until IBD_Pool's destructor runs.
The malloc'd slabs just grow as needed and are freed at the end.


### Why the stack reversal in getUnprocessed?

This confused me at first. The stack stores nodes in reverse genomic
order (newest position on top). But when we need to reprocess the
tail nodes, we need to feed them back into add_node in genomic order
(oldest first). That's because add_node pushes onto the stack, and
the stack needs to maintain its "newest on top" invariant.

The reverse() call flips the stack to genomic order, then the split
at `end` gives us the tail in the right order to pop-and-reprocess.

```
NORMAL STACK ORDER (reverse genomic):
  top -> [newest] -> ... -> [end] -> ... -> [oldest/start]

AFTER reverse() (genomic order):
  top -> [oldest] -> ... -> [end] -> ... -> [newest/start]

SPLIT at end:
  segment:     [oldest] -> ... -> [end]     <- gets reclaimed
  unprocessed: [end+1]  -> ... -> [newest]  <- popped one by one into add_node
```


### Worked example: getUnprocessed producing a second segment

It's hard to construct this case intuitively because the cumulative LOD
has to go negative (triggering reachedMax) while the tail after `end`
still contains enough positive LOD to form a new segment above threshold.
The key is a strong first peak, a big dip that doesn't quite kill the
cumulative, positive sites that never catch up to the first peak's
cumulative, and then a crash.

Sites at positions 100-600, LODs: +2, +8, -8, +4, +3, -15, threshold=3.

#### Step 1: Push pos=100, lod=+2

```
segment: [100: lod=+2, cum=2]
          ^top ^start ^end
```
Singleton. Initialize recorders.

#### Step 2: Push pos=200, lod=+8

```
segment: [200: lod=+8, cum=10] -> [100: cum=2]
          ^top                      ^start
```
topIsNewMax (10 >= 2) -> setEnd=200, reclaim interior (none, adjacent).
```
segment: [200: cum=10] -> [100]
          ^top ^end         ^start
```

#### Step 3: Push pos=300, lod=-8

```
segment: [300: lod=-8, cum=2] -> [200: cum=10] -> [100]
          ^top                     ^end              ^start
```
topIsNewMax? 2 >= 10 -> NO.  reachedMax? 2 < 0 -> NO.
Node 300 just sits on top, dragging the sum down but not enough.

#### Step 4: Push pos=400, lod=+4

```
segment: [400: lod=+4, cum=6] -> [300: cum=2] -> [200: cum=10] -> [100]
          ^top                                      ^end             ^start
```
topIsNewMax? 6 >= 10 -> NO. The +4 helped but can't overcome the -8 dip.
reachedMax? 6 < 0 -> NO.

#### Step 5: Push pos=500, lod=+3

```
segment: [500: lod=+3, cum=9] -> [400: cum=6] -> [300: cum=2] -> [200: cum=10] -> [100]
          ^top                                                      ^end             ^start
```
topIsNewMax? 9 >= 10 -> NO. Still can't catch the first peak!
reachedMax? 9 < 0 -> NO.

Note: positions 400 and 500 have strong positive LODs (+4, +3) but their
cumulative from start (6, 9) never reaches end's cumulative (10). They're
"trapped" in the tail of a dominant first peak.

#### Step 6: Push pos=600, lod=-15

```
segment: [600: lod=-15, cum=-6] -> [500: cum=9] -> [400: cum=6] -> [300: cum=2] -> [200: cum=10] -> [100]
          ^top                                                        ^end             ^start
```
topIsNewMax? -6 >= 10 -> NO.
reachedMax? -6 < 0 -> YES!

**6a. Emit segment 1:**
  endLod=10 >= threshold(3) -> YES!
  Output: `sample  chrom  100  200  10.0`

**6b. getUnprocessed():**

Reverse the stack to genomic order:
```
BEFORE: [600] -> [500] -> [400] -> [300] -> [200] -> [100]
                                              ^end     ^start
AFTER:  [100] -> [200] -> [300] -> [400] -> [500] -> [600]
                  ^end                                 ^start
```

Split at end:
```
segment:     [100] -> [200]     (start=end=[200])
unprocessed: [300] -> [400] -> [500] -> [600]
                                         ^start
```

**6c. reclaim_stack(&segment):** nodes [100] and [200] go to pool.

**6d. Reprocess unprocessed tail:**

Pop [300] (lod=-8): segment is empty AND lod < 0 -> reclaim, skip.

Pop [400] (lod=+4): segment is empty, lod >= 0 -> PUSH!
```
segment: [400: lod=+4, cum=4]
          ^top ^start ^end
```
Singleton, fresh segment started from the tail!

Pop [500] (lod=+3): push.
```
segment: [500: lod=+3, cum=7] -> [400: cum=4]
          ^top                     ^start
```
topIsNewMax? 7 >= 4 -> YES! setEnd=500.
```
segment: [500: cum=7] -> [400]
          ^top ^end        ^start
```

Pop [600] (lod=-15): push.
```
segment: [600: lod=-15, cum=-8] -> [500: cum=7] -> [400]
          ^top                       ^end             ^start
```
topIsNewMax? -8 >= 7 -> NO.
reachedMax? -8 < 0 -> YES!

**Emit segment 2:**
  endLod=7 >= threshold(3) -> YES!
  Output: `sample  chrom  400  500  7.0`

**Final result: TWO segments emitted from the same stretch of data:**
```
Segment 1: pos 100-200, LOD=10  (from the original stack)
Segment 2: pos 400-500, LOD=7   (from reprocessing the tail!)
```

Without getUnprocessed, the nodes at positions 400 and 500 would have
been thrown away with the rest of the stack when segment 1 was emitted.
Those +4 and +3 LOD scores would have been lost, and the 400-500
introgressed region would never have been reported.

This happens biologically when two introgressed tracts are close together
with a short non-introgressed gap between them. The first tract's LOD
peak dominates, but after it's emitted, the second tract's LODs are
strong enough to stand on their own.


### Summary of what I didn't understand before

1. **IBD_Pool** is a free-list memory allocator. It pre-allocates IBD_Node
   objects in bulk and recycles them. All samples share one pool. The pool's
   internal storage is itself an IBD_Stack (just used as a simple linked list).

2. **Introgressed regions are found via max-subarray**: accumulate LOD scores
   as a running sum, track the maximum, emit when the sum goes negative (if
   max exceeded threshold), then reprocess the tail for potential adjacent
   segments. The three reclaim methods handle recycling nodes at different
   stages of this process.