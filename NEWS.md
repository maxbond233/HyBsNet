# HyBsNet 0.2.0.9000

* Adds topology-aware import, freezing, comparison, and reconciliation of
  publication layouts. Shared coordinates can remain exact while replacement
  nodes are explicitly anchored to positioned neighbours.
* Adds standardized differential-gene, expression, pathway, and custom node
  overlays with separate tested, available, and significant states. Missing
  values are never silently converted to zero.
* Adds fixed-layout overlay plotting for differential expression, GSEA NES,
  pathway-activity scores, and other node-level metrics.
* Adds reusable ego, explicit-node, pattern, community, and group subnetwork
  queries, with explicit induced versus seed-incident edge semantics.
* Adds frozen-parent and compact subnetwork layouts, local plots, context plots,
  a Source-able example workflow, and real-data regression validation.
* The Source-able example now accepts explicit DEG, edge, frozen-layout,
  optional pathway, and output paths instead of inferring a HyBs project
  directory structure.

# HyBsNet 0.1.0.9000

* Initial internal package for the HyBs Figure 3 network workflow.
* Separates DEG/Jaccard/network analysis, deterministic layout, figure styling,
  and file export.
* Adds current Figure 3 and colour-blind-friendly style presets.
* Adds unit tests and a real-data regression script for Figure 3A/B/F/G.
