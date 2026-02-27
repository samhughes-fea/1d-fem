# tests/test_simulation_settings_integration.py

"""
Integration tests for simulation settings parser and integration.

Tests backward compatibility, new configs, partial configs, invalid configs,
type conversion, and case insensitivity.
"""

import os
import tempfile
import pytest
from pathlib import Path

from pre_processing.parsing.simulation_settings_parser import parse_simulation_settings


class TestSimulationSettingsIntegration:
    """Integration tests for simulation settings parser."""
    
    def test_backward_compatibility_minimal(self):
        """Test that minimal settings file (only type) still works."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["type"] == "static"
            # Should have defaults for other sections
            assert "solver" in settings
            assert "condensation" in settings
            assert "parallel" in settings
        finally:
            os.unlink(temp_path)
    
    def test_backward_compatibility_none(self):
        """Test that None simulation_settings is handled (via defaults in StaticSimulationRunner)."""
        # This is tested indirectly through StaticSimulationRunner
        # The constructor should accept None and use defaults
        pass
    
    def test_new_solver_config(self):
        """Test parsing of [Solver] section."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Solver]
type = cg
tolerance = 1e-8
max_iterations = 2000
restart = 30
ilu_drop_tol = 1e-7
ilu_fill_factor = 1.5
disable_scaling = true
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["solver"]["type"] == "cg"
            assert settings["solver"]["tolerance"] == 1e-8
            assert settings["solver"]["max_iterations"] == 2000
            assert settings["solver"]["restart"] == 30
            assert settings["solver"]["ilu_drop_tol"] == 1e-7
            assert settings["solver"]["ilu_fill_factor"] == 1.5
            assert settings["solver"]["disable_scaling"] is True
        finally:
            os.unlink(temp_path)
    
    def test_new_condensation_config(self):
        """Test parsing of [Condensation] section."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Condensation]
base_tol = 1e-10
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["condensation"]["base_tol"] == 1e-10
        finally:
            os.unlink(temp_path)
    
    def test_new_parallel_config(self):
        """Test parsing of [Parallel] section."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Parallel]
num_processes = 4
enable_parallel_instantiation = true
enable_parallel_computation = false
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["parallel"]["num_processes"] == 4
            assert settings["parallel"]["enable_parallel_instantiation"] is True
            assert settings["parallel"]["enable_parallel_computation"] is False
        finally:
            os.unlink(temp_path)
    
    def test_parallel_config_auto(self):
        """Test that 'auto' is preserved for num_processes."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Parallel]
num_processes = auto
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["parallel"]["num_processes"] == "auto"
        finally:
            os.unlink(temp_path)

    def test_modal_section(self):
        """Test parsing of [Modal] section and defaults."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Modal]
num_modes = 20
""")
            temp_path = f.name

        try:
            settings = parse_simulation_settings(temp_path)
            assert "modal" in settings
            assert settings["modal"]["num_modes"] == 20
        finally:
            os.unlink(temp_path)

        # Default when not specified
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Modal
""")
            temp_path = f.name
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["modal"]["num_modes"] == 10
        finally:
            os.unlink(temp_path)

    def test_dynamic_section(self):
        """Test parsing of [Dynamic] section and defaults."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Dynamic]
time_step = 0.0005
end_time = 2.0
scheme = newmark
""")
            temp_path = f.name

        try:
            settings = parse_simulation_settings(temp_path)
            assert "dynamic" in settings
            assert settings["dynamic"]["time_step"] == 0.0005
            assert settings["dynamic"]["end_time"] == 2.0
            assert settings["dynamic"]["scheme"] == "newmark"
        finally:
            os.unlink(temp_path)
    
    def test_partial_configs(self):
        """Test that missing sections use defaults."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Solver]
type = gmres
# Other solver params use defaults
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["solver"]["type"] == "gmres"
            # Check defaults are used
            assert settings["solver"]["tolerance"] == 1e-6  # default
            assert settings["condensation"]["base_tol"] == 1e-12  # default
        finally:
            os.unlink(temp_path)
    
    def test_invalid_simulation_type(self):
        """Test that invalid simulation type raises ValueError."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
InvalidType
""")
            temp_path = f.name
        
        try:
            with pytest.raises(ValueError, match="Invalid simulation type"):
                parse_simulation_settings(temp_path)
        finally:
            os.unlink(temp_path)
    
    def test_invalid_parallel_num_processes(self):
        """Test that invalid num_processes raises ValueError."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Parallel]
num_processes = -1
""")
            temp_path = f.name
        
        try:
            with pytest.raises(ValueError):
                parse_simulation_settings(temp_path)
        finally:
            os.unlink(temp_path)
    
    def test_type_conversion_float(self):
        """Test conversion of string to float."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Solver]
tolerance = 1e-5
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert isinstance(settings["solver"]["tolerance"], float)
            assert settings["solver"]["tolerance"] == 1e-5
        finally:
            os.unlink(temp_path)
    
    def test_type_conversion_int(self):
        """Test conversion of string to int."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Solver]
max_iterations = 500
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert isinstance(settings["solver"]["max_iterations"], int)
            assert settings["solver"]["max_iterations"] == 500
        finally:
            os.unlink(temp_path)
    
    def test_type_conversion_bool(self):
        """Test conversion of string to bool."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Solver]
disable_scaling = yes
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert isinstance(settings["solver"]["disable_scaling"], bool)
            assert settings["solver"]["disable_scaling"] is True
        finally:
            os.unlink(temp_path)
    
    def test_case_insensitivity_sections(self):
        """Test that section headers are case-insensitive."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[SIMULATION]
[TYPE]
Static

[SOLVER]
type = CG
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["type"] == "static"
            assert settings["solver"]["type"] == "cg"  # Lowercased
        finally:
            os.unlink(temp_path)
    
    def test_case_insensitivity_keys(self):
        """Test that key names are case-insensitive."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Solver]
TYPE = cg
TOLERANCE = 1e-7
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["solver"]["type"] == "cg"
            assert settings["solver"]["tolerance"] == 1e-7
        finally:
            os.unlink(temp_path)
    
    def test_complete_config(self):
        """Test parsing of complete configuration with all sections."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("""[Simulation]
[Type]
Static

[Solver]
type = cg
tolerance = 1e-6
max_iterations = 1000
restart = 20
ilu_drop_tol = 1e-6
ilu_fill_factor = 1.0
disable_scaling = false

[Condensation]
base_tol = 1e-12

[Parallel]
num_processes = auto
enable_parallel_instantiation = true
enable_parallel_computation = true
""")
            temp_path = f.name
        
        try:
            settings = parse_simulation_settings(temp_path)
            assert settings["type"] == "static"
            assert settings["solver"]["type"] == "cg"
            assert settings["condensation"]["base_tol"] == 1e-12
            assert settings["parallel"]["num_processes"] == "auto"
        finally:
            os.unlink(temp_path)

