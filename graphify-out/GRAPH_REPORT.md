# Graph Report - .  (2026-08-10)

## Corpus Check
- Large corpus: 141 files · ~895,914 words. Semantic extraction will be expensive (many Claude tokens). Consider running on a subfolder.

## Summary
- 495 nodes · 771 edges · 45 communities (35 shown, 10 thin omitted)
- Extraction: 99% EXTRACTED · 1% INFERRED · 0% AMBIGUOUS · INFERRED: 11 edges (avg confidence: 0.81)
- Token cost: 55,663 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Solid FEM Solver Core|Solid FEM Solver Core]]
- [[_COMMUNITY_Assembly & Projection Routines|Assembly & Projection Routines]]
- [[_COMMUNITY_Wolf Git Hooks (JS)|Wolf Git Hooks (JS)]]
- [[_COMMUNITY_Wolf Config|Wolf Config]]
- [[_COMMUNITY_Per-Task Data Buffers|Per-Task Data Buffers]]
- [[_COMMUNITY_Density & Material Models|Density & Material Models]]
- [[_COMMUNITY_Build & Dependencies|Build & Dependencies]]
- [[_COMMUNITY_GUI Info Dialogs|GUI Info Dialogs]]
- [[_COMMUNITY_Scratch Data  FE Values|Scratch Data / FE Values]]
- [[_COMMUNITY_OpenWolf Protocol Docs|OpenWolf Protocol Docs]]
- [[_COMMUNITY_Visualization Video Generators|Visualization Video Generators]]
- [[_COMMUNITY_UQPH Scratch Fields|UQPH Scratch Fields]]
- [[_COMMUNITY_ScratchData_K Tensors|ScratchData_K Tensors]]
- [[_COMMUNITY_ScratchData_RHS Face Values|ScratchData_RHS Face Values]]
- [[_COMMUNITY_3D Rotation Utility|3D Rotation Utility]]
- [[_COMMUNITY_Claude Hooks Settings|Claude Hooks Settings]]
- [[_COMMUNITY_Cron Engine State|Cron Engine State]]
- [[_COMMUNITY_Initial Density Condition|Initial Density Condition]]
- [[_COMMUNITY_Design QC Report|Design QC Report]]
- [[_COMMUNITY_GUI Default Values|GUI Default Values]]
- [[_COMMUNITY_SQPH Scratch Data|SQPH Scratch Data]]
- [[_COMMUNITY_Parameter Header|Parameter Header]]
- [[_COMMUNITY_Bug Log|Bug Log]]
- [[_COMMUNITY_Cron Manifest|Cron Manifest]]
- [[_COMMUNITY_Suggestions Store|Suggestions Store]]
- [[_COMMUNITY_SQPH Global Copy|SQPH Global Copy]]
- [[_COMMUNITY_UQPH Global Copy|UQPH Global Copy]]
- [[_COMMUNITY_Dylib Fix Script|Dylib Fix Script]]
- [[_COMMUNITY_NPM Package Manifest|NPM Package Manifest]]
- [[_COMMUNITY_Assistant Identity Config|Assistant Identity Config]]

## God Nodes (most connected - your core abstractions)
1. `Solid` - 99 edges
2. `run()` - 17 edges
3. `Solid<dim>::ScratchData_UQPH` - 17 edges
4. `BlockVector` - 16 edges
5. `solve_nonlinear_timestep()` - 16 edges
6. `main()` - 15 edges
7. `Solid<dim>::PerTaskData_SC` - 15 edges
8. `Solid<dim>::ScratchData_RHS` - 15 edges
9. `getWolfDir()` - 14 edges
10. `ensureWolfDir()` - 14 edges

## Surprising Connections (you probably didn't know these)
- `OpenWolf Protocol Enforcement Rules` --semantically_similar_to--> `OpenWolf Operating Protocol`  [INFERRED] [semantically similar]
  .claude/rules/openwolf.md → .wolf/OPENWOLF.md
- `Python3 User Interface Component` --references--> `customtkinter`  [INFERRED]
  README.md → requirements.txt
- `Python3 User Interface Component` --references--> `numpy`  [INFERRED]
  README.md → requirements.txt
- `Python3 User Interface Component` --references--> `opencv-python`  [INFERRED]
  README.md → requirements.txt
- `Python3 User Interface Component` --references--> `py2app`  [INFERRED]
  README.md → requirements.txt

## Import Cycles
- None detected.

## Hyperedges (group relationships)
- **OpenWolf Cross-Session Context Management** — wolf_openwolf_protocol, wolf_anatomy_map, wolf_cerebrum_memory, wolf_memory_log, wolf_openwolf_buglog [EXTRACTED 0.90]
- **BFSimulator Two-Component Architecture** — readme_bfsimulator, readme_dealii, readme_python_interface, cmakelists_brain_growth [EXTRACTED 0.85]

## Communities (45 total, 10 thin omitted)

### Community 0 - "Solid FEM Solver Core"
Cohesion: 0.02
Nodes (87): AffineConstraints, BlockSparseMatrix, BlockSparsityPattern, determine_component_extractors(), Solid, assemble_sc, assemble_sc_one_cell, assemble_system_rhs (+79 more)

### Community 1 - "Assembly & Projection Routines"
Cohesion: 0.07
Nodes (57): BlockVector, assemble_sc(), assemble_sc_one_cell(), assemble_system_rhs(), assemble_system_rhs_one_cell(), assemble_system_tangent(), assemble_system_tangent_one_cell(), PointHistory (+49 more)

### Community 3 - "Wolf Git Hooks (JS)"
Cohesion: 0.18
Nodes (34): main(), autoDetectBugFix(), detectFixPattern(), extractCalls(), extractChangedLines(), extractCSSProps(), findOperatorChange(), main() (+26 more)

### Community 4 - "Wolf Config"
Cohesion: 0.05
Nodes (39): auto_scan_on_init, exclude_patterns, max_description_length, max_files, rescan_interval_hours, max_tokens, reflection_frequency, api_key_env (+31 more)

### Community 5 - "Per-Task Data Buffers"
Cohesion: 0.10
Nodes (17): Solid<dim>::PerTaskData_K, cell_matrix, local_dof_indices, Solid<dim>::PerTaskData_RHS, local_dof_indices, Solid<dim>::PerTaskData_SC, k_B, k_bar (+9 more)

### Community 6 - "Density & Material Models"
Cohesion: 0.16
Nodes (15): compute_denisty_source(), compute_denisty_source_ORGC(), get_dcc_r(), get_v_r(), Point, heaviside_function(), namespace, NonStandardTensors() (+7 more)

### Community 7 - "Build & Dependencies"
Cohesion: 0.16
Nodes (14): Brain_growth CMake Target, deal.II Library Dependency (CMake), fix_dylib_deps.sh Post-Build Command, BFSimulator Project, deal.II C++ Computational Model, eLife Cortical Folding OSVZ Paper, Paraview Visualization Tool, Python3 User Interface Component (+6 more)

### Community 8 - "GUI Info Dialogs"
Cohesion: 0.18
Nodes (11): About_author(), About_programm(), c_k_info(), callback(), Copy_right(), degree_info(), nonlinear_it_info(), poisson_ratio_info() (+3 more)

### Community 9 - "Scratch Data / FE Values"
Cohesion: 0.29
Nodes (7): FiniteElement, QGauss, ScratchData_K, ScratchData_RHS, ScratchData_SQPH, ScratchData_UQPH, UpdateFlags

### Community 10 - "OpenWolf Protocol Docs"
Cohesion: 0.35
Nodes (11): CLAUDE.md OpenWolf Bootstrap, OpenWolf Protocol Enforcement Rules, Project Anatomy File Map, Do-Not-Repeat List, Cerebrum Learning Memory, Chronological Action Memory Log, Bug Logging Mechanism (buglog.json), Design QC Screenshot Workflow (+3 more)

### Community 12 - "UQPH Scratch Fields"
Cohesion: 0.20
Nodes (10): Solid<dim>::ScratchData_UQPH, fe_values_ref, old_old_solution_grads_u, old_old_solution_value_c, old_solution_grads_u, old_solution_value_c, solution_grads_c, solution_grads_u (+2 more)

### Community 13 - "ScratchData_K Tensors"
Cohesion: 0.22
Nodes (8): Tensor, Solid<dim>::ScratchData_K, fe_values_ref, grad_Nx, grad_Nx_c, Nx, symm_grad_Nx, SymmetricTensor

### Community 14 - "ScratchData_RHS Face Values"
Cohesion: 0.22
Nodes (8): Solid<dim>::ScratchData_RHS, fe_face_values_ref, fe_values_ref, grad_Nx, grad_Nx_c, Nx, symm_grad_Nx, FEFaceValues

### Community 15 - "3D Rotation Utility"
Cohesion: 0.29
Nodes (5): dim, Point, Rotate3d, angle, axis

### Community 16 - "Claude Hooks Settings"
Cohesion: 0.33
Nodes (5): hooks, PostToolUse, PreToolUse, SessionStart, Stop

### Community 17 - "Cron Engine State"
Cohesion: 0.33
Nodes (5): dead_letter_queue, engine_status, execution_log, last_heartbeat, upcoming

### Community 18 - "Initial Density Condition"
Cohesion: 0.40
Nodes (4): InitialValueC, dvision_raduis, dvision_value, Function

### Community 19 - "Design QC Report"
Cohesion: 0.40
Nodes (4): captured_at, captures, estimated_tokens, total_size_kb

### Community 20 - "GUI Default Values"
Cohesion: 0.50
Nodes (4): messageWindow(), set_2D_default(), set_3D_default(), set_default_values()

### Community 21 - "SQPH Scratch Data"
Cohesion: 0.50
Nodes (3): Solid<dim>::ScratchData_SQPH, fe_values_ref, FEValues

## Knowledge Gaps
- **191 isolated node(s):** `SessionStart`, `PreToolUse`, `PostToolUse`, `Stop`, `version` (+186 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **10 thin communities (<3 nodes) omitted from report** — run `graphify query` to explore isolated nodes.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `Solid` connect `Solid FEM Solver Core` to `Assembly & Projection Routines`, `Per-Task Data Buffers`, `Scratch Data / FE Values`, `3D Rotation Utility`, `SQPH Global Copy`, `UQPH Global Copy`?**
  _High betweenness centrality (0.143) - this node is a cross-community bridge._
- **Why does `Solid<dim>::PerTaskData_SC` connect `Per-Task Data Buffers` to `Assembly & Projection Routines`?**
  _High betweenness centrality (0.020) - this node is a cross-community bridge._
- **Why does `Solid<dim>::ScratchData_UQPH` connect `UQPH Scratch Fields` to `Assembly & Projection Routines`, `Per-Task Data Buffers`, `Scratch Data / FE Values`, `ScratchData_K Tensors`, `3D Rotation Utility`, `SQPH Scratch Data`?**
  _High betweenness centrality (0.019) - this node is a cross-community bridge._
- **What connects `SessionStart`, `PreToolUse`, `PostToolUse` to the rest of the system?**
  _191 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Solid FEM Solver Core` be split into smaller, more focused modules?**
  _Cohesion score 0.022988505747126436 - nodes in this community are weakly interconnected._
- **Should `Assembly & Projection Routines` be split into smaller, more focused modules?**
  _Cohesion score 0.07075873827791987 - nodes in this community are weakly interconnected._
- **Should `Python GUI Validation` be split into smaller, more focused modules?**
  _Cohesion score 0.03773584905660377 - nodes in this community are weakly interconnected._