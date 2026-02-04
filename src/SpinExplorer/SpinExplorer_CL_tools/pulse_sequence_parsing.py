from pathlib import Path
from typing import Optional, Dict, List
from experiment_config import ExperimentConfigStore


class PulseSequenceError(Exception):
    """Base exception for pulse sequence parsing errors."""
    pass


class PulseProgramNotFoundError(PulseSequenceError):
    """Raised when pulseprogram file is not found."""
    pass


class NoConfigurationError(PulseSequenceError):
    """Raised when no configuration is found for a pulse sequence."""
    pass


class ConfigurationRegistry:
    """
    Registry linking pulse sequence names to their processing configurations.
    Supports one-to-many relationships (one sequence can have multiple configs).
    """
    
    def __init__(self):
        self._registry: Dict[str, List[ExperimentConfigStore]] = {}
    
    def register(self, sequence_name: str, config: ExperimentConfigStore) -> None:
        """
        Register a configuration for a pulse sequence.
        
        Args:
            sequence_name: Name of the pulse sequence (e.g., 'hsqcetfpf3gp')
            config: ProcessingConfig object
        """
        if sequence_name not in self._registry:
            self._registry[sequence_name] = []
        self._registry[sequence_name].append(config)
    
    def get_configs(self, sequence_name: str) -> List[ExperimentConfigStore]:
        """
        Get all configurations for a pulse sequence.
        
        Args:
            sequence_name: Name of the pulse sequence
            
        Returns:
            List of ProcessingConfig objects
            
        Raises:
            NoConfigurationError: If no configuration exists for this sequence
        """
        if sequence_name not in self._registry:
            raise NoConfigurationError(
                f"No configuration found for pulse sequence '{sequence_name}'. "
                f"Available sequences: {list(self._registry.keys())}"
            )
        return self._registry[sequence_name]
    
    def get_default_config(self, sequence_name: str) -> ExperimentConfigStore:
        """
        Get the first (default) configuration for a pulse sequence.
        
        Args:
            sequence_name: Name of the pulse sequence
            
        Returns:
            ProcessingConfig object
            
        Raises:
            NoConfigurationError: If no configuration exists for this sequence
        """
        configs = self.get_configs(sequence_name)
        return configs[0]
    
    def has_config(self, sequence_name: str) -> bool:
        """Check if a configuration exists for a pulse sequence."""
        return sequence_name in self._registry
    
    def list_sequences(self) -> List[str]:
        """List all registered pulse sequences."""
        return list(self._registry.keys())


class PulseSequenceParser:
    """
    Parser for Bruker pulse sequence files.
    
    Extracts pulse sequence name from 'pulseprogram' files.
    """
    
    def __init__(self, folder: Optional[Path] = None, 
                 config_registry: Optional[ConfigurationRegistry] = None):
        """
        Initialize the parser.
        
        Args:
            folder: Path to folder containing 'pulseprogram' file.
                   Defaults to current working directory.
            config_registry: Optional ConfigurationRegistry for validation.
                           If provided, will check that parsed sequences have configs.
        """
        self.folder = Path(folder) if folder else Path.cwd()
        self.config_registry = config_registry
        self._sequence_name: Optional[str] = None
    
    def parse(self) -> str:
        """
        Parse the pulseprogram file and extract the sequence name.
        
        Returns:
            The pulse sequence name (e.g., 'hsqcetfpf3gp')
            
        Raises:
            PulseProgramNotFoundError: If 'pulseprogram' file doesn't exist
            PulseSequenceError: If no valid sequence name found
            NoConfigurationError: If config_registry provided and no config exists
        """
        pulseprogram_path = self.folder / 'pulseprogram'
        
        if not pulseprogram_path.exists():
            raise PulseProgramNotFoundError(
                f"'pulseprogram' file not found in {self.folder}"
            )
        
        # Read file and find first line starting with semicolon
        with open(pulseprogram_path, 'r', encoding='utf-8') as f:
            for line in f:
                stripped = line.strip()
                if stripped.startswith(';'):
                    # Extract everything after the semicolon
                    sequence_name = stripped[1:].strip()
                    if sequence_name:  # Make sure it's not empty
                        self._sequence_name = sequence_name
                        
                        # Validate against config registry if provided
                        if self.config_registry is not None:
                            if not self.config_registry.has_config(sequence_name):
                                raise NoConfigurationError(
                                    f"No configuration found for pulse sequence '{sequence_name}'"
                                )
                        
                        return sequence_name
        
        raise PulseSequenceError(
            f"No valid pulse sequence name found in {pulseprogram_path}"
        )
    
    @property
    def sequence_name(self) -> str:
        """
        Get the parsed sequence name.
        
        Returns:
            The sequence name
            
        Raises:
            PulseSequenceError: If parse() hasn't been called yet
        """
        if self._sequence_name is None:
            raise PulseSequenceError("No sequence parsed yet. Call parse() first.")
        return self._sequence_name
    
    def get_config(self) -> ExperimentConfigStore:
        """
        Get the default configuration for the parsed sequence.
        
        Returns:
            ProcessingConfig object
            
        Raises:
            PulseSequenceError: If parse() hasn't been called or no registry provided
            NoConfigurationError: If no configuration exists
        """
        if self._sequence_name is None:
            raise PulseSequenceError("No sequence parsed yet. Call parse() first.")
        
        if self.config_registry is None:
            raise PulseSequenceError("No configuration registry provided to parser.")
        
        return self.config_registry.get_default_config(self._sequence_name)
    
    def get_all_configs(self) -> List[ExperimentConfigStore]:
        """
        Get all configurations for the parsed sequence.
        
        Returns:
            List of ProcessingConfig objects
            
        Raises:
            PulseSequenceError: If parse() hasn't been called or no registry provided
            NoConfigurationError: If no configuration exists
        """
        if self._sequence_name is None:
            raise PulseSequenceError("No sequence parsed yet. Call parse() first.")
        
        if self.config_registry is None:
            raise PulseSequenceError("No configuration registry provided to parser.")
        
        return self.config_registry.get_configs(self._sequence_name)


# Example usage:
if __name__ == "__main__":
    # Create a configuration registry
    registry = ConfigurationRegistry()
    
    # Register some configurations
    hsqc_config = ProcessingConfig(
        name="hsqc_standard",
        description="Standard HSQC processing"
    )
    registry.register("hsqcetfpf3gp", hsqc_config)
    
    # You can register multiple configs for the same sequence
    hsqc_config_alt = ProcessingConfig(
        name="hsqc_highres",
        description="High resolution HSQC processing"
    )
    registry.register("hsqcetfpf3gp", hsqc_config_alt)
    
    # Parse a pulse sequence
    try:
        parser = PulseSequenceParser(
            folder=Path("/path/to/your/data"),
            config_registry=registry
        )
        sequence_name = parser.parse()
        print(f"Found pulse sequence: {sequence_name}")
        
        # Get the default config
        config = parser.get_config()
        print(f"Using config: {config.name}")
        
        # Or get all configs
        all_configs = parser.get_all_configs()
        print(f"Available configs: {[c.name for c in all_configs]}")
        
    except PulseSequenceError as e:
        print(f"Error: {e}")