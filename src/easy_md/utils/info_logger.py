"""
Logging configuration utility for the application.

This module provides centralized logging configuration and logger instance
creation for the entire application. It sets up file-based logging with
appropriate formatting and level configuration.
"""

import logging
from logging.handlers import RotatingFileHandler
import os
from pathlib import Path

# Define log file path
current_file = Path(__file__)
project_root = current_file.parent.parent.parent.parent  # Go up from utils/easy_md/src/easy-md/
logs_dir = project_root / 'logs'
log_file_path = logs_dir / 'info_eror.log'

# Create logs directory if it doesn't exist
os.makedirs(logs_dir, exist_ok=True)

# Create rotating file handler
rotating_handler = RotatingFileHandler(
    filename=log_file_path,
    maxBytes=5_000_000,      # Rotate after 5 MB
    backupCount=3            # Keep up to 3 old log files
)

formatter = logging.Formatter(
    '%(asctime)s - %(levelname)s - [%(name)s:%(lineno)d] - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
rotating_handler.setFormatter(formatter)

# Create console handler
console_handler = logging.StreamHandler()
console_handler.setFormatter(formatter)

# Apply logging configuration to all loggers
logging.basicConfig(level=logging.INFO, handlers=[rotating_handler, console_handler])