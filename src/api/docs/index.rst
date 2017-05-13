===========================
Gromacs External Interfaces
===========================

..  toctree::
    :maxdepth: 2
    :caption: Documentation sections

    intro
    userguide
    developerguide
    components
    layers
    reference
    changelog

..
  Documentation walk-through:

  * Description: goals and scope of the project.
  * User guides: important concepts, how to get started, "How-to"
    links to documentation by example begin as user stories.
    Cross-link with components doc (feature) when possible.
  * Developer guides: further important concepts, design principles,
    architecture-level requirements, protocols and interfaces
    (sequence diagrams, object diagrams, class diagrams (where
    available))
  * Components: major sections of the project, hierarchical
    functionality describes features to guide development.
  * Layers: interface / API layers. Use cases and scenarios
    expressed here help identify features and clarify implementation levels
    for the various components.
  * Reference: auto-generated documentation for high level interfaces
    and lower level API for Python and C++.

  Terminology conventions:

  "User stories" are very concise expression of features from a user's perspective, as favored in Agile development. We can use this format as a uniform way to transcribe feature requests and, in the absence of a Scrum system, express them in the User Guide. When implemented, the user story becomes a link to a "how to" or other relevant documentation.

  "Use case" is a set of scenarios describing a user goal, used to capture functional requirements. Use cases inform feature targets for development iterations. Use cases are expressed at a specific level of interaction or API layer, but a step in a use case can be another use case.

  "Scenario" is a sequence of steps for a client to carry out a task.
  Exceptions describe alternate paths.

  "Feature" expresses the target of a development iteration targeting use case(s) from the developer perspective.

  "Requirements" are clear and verifiable design constraints.
  Functional requirements are implementation-agnostic specifications on a feature, expressed by the developer as interpreted from a targeted use case, that can be verified to show the API does what it claims.
