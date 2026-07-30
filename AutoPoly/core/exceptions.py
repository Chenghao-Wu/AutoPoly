"""
Custom exception classes for AutoPoly.

This module defines the exception hierarchy used throughout AutoPoly for proper
error handling and testing.

Exceptions are used instead of sys.exit() to enable:
- Proper unit testing with pytest.raises()
- Graceful cleanup and error recovery
- Full stack traces for debugging
- Programmatic error handling by library users
"""


class AutoPolyError(Exception):
    """Base exception for all AutoPoly errors.

    All AutoPoly-specific exceptions inherit from this class, allowing
    users to catch all AutoPoly errors with a single except clause.

    Example:
        try:
            poly = Polymer(...)
        except AutoPolyError as e:
            print(f"AutoPoly error: {e}")
    """
    pass


class ValidationError(AutoPolyError):
    """Raised when input validation fails.

    This exception is raised when user-provided input does not meet
    validation requirements, such as:
    - Invalid parameter values (e.g., ChainNum=0)
    - Invalid sequence length
    - Missing required fields
    - Invalid SMILES strings

    Example:
        >>> from AutoPoly.core.exceptions import ValidationError
        >>> try:
        ...     Polymer(ChainNum=0, Sequence=["[*]CC[*]"])
        ... except ValidationError as e:
        ...     print(f"Validation failed: {e}")
    """

    pass


class GenerationError(AutoPolyError):
    """Raised when monomer or polymer generation fails.

    This exception is raised when the generation process encounters
    an error, such as:
    - RDKit monomer generation failures
    - SMILES parsing errors
    - Force field application failures
    - Moltemplate file generation errors

    Example:
        >>> from AutoPoly.core.exceptions import GenerationError
        >>> try:
        ...     generate_monomer_from_smiles("invalid_smiles")
        ... except GenerationError as e:
        ...     print(f"Generation failed: {e}")
    """

    pass


class WorkflowError(AutoPolyError):
    """Raised when workflow execution fails.

    This exception is raised when the overall workflow process fails,
    such as:
    - Moltemplate execution failures
    - LAMMPS data file generation errors
    - File system operation failures
    - External tool integration errors

    Example:
        >>> from AutoPoly.core.exceptions import WorkflowError
        >>> try:
        ...     generate(system, name, models, force_field="oplsaa")
        ... except WorkflowError as e:
        ...     print(f"Workflow failed: {e}")
    """

    pass
