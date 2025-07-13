"""
Logging configuration utility for the application.

This module provides centralized logging configuration and logger instance
creation for the entire application. It sets up file-based logging with
appropriate formatting and level configuration.
"""

import logging
from logging.handlers import RotatingFileHandler
import os

# Define log file path
app_dir = os.path.dirname(os.path.dirname(__file__))  # Go up from utils/ to easy_md/
log_file_path = os.path.join(app_dir, 'info_eror.log')

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