import logging
from logging.handlers import RotatingFileHandler
import os
import sys
from pathlib import Path
from datetime import datetime
from typing import Optional, Any


class LoggerHelper:
    _configured: bool = False
    _logger: Optional[logging.Logger] = None

    @classmethod
    def configure(
        cls,
        name: str = "trial-statis-server",
        level: Optional[str] = None,
        log_dir: str = "logs",
        max_bytes: int = 10 * 1024 * 1024,
        backup_count: int = 5,
        fmt: Optional[str] = None,
    ) -> None:
        """Configure the logger to write logs directly in logs directory with yyyy-MM-dd.log naming.
        
        Args:
            name: Logger name
            level: Log level (defaults to INFO)
            log_dir: Base directory for logs (default: "logs")
            max_bytes: Maximum size per log file in bytes
            backup_count: Number of backup files to keep
            fmt: Log format string
        """
        if cls._configured and cls._logger is not None:
            return

        # Set log level
        level = level or os.getenv("LOG_LEVEL", "INFO")
        numeric = getattr(logging, level.upper(), logging.INFO)

        # Create logger
        logger = logging.getLogger(name)
        logger.setLevel(numeric)

        # If handlers already present, assume configured externally
        if logger.handlers:
            cls._logger = logger
            cls._configured = True
            return

        # Create logs directory if it doesn't exist
        logs_path = Path(log_dir)
        logs_path.mkdir(parents=True, exist_ok=True)

        # Generate filename with yyyy-MM-dd format
        today = datetime.now().strftime("%Y-%m-%d")
        log_filename = f"{today}.log"
        file_path = logs_path / log_filename

        # Create file handler with rotation
        file_handler = RotatingFileHandler(
            str(file_path), 
            maxBytes=max_bytes, 
            backupCount=backup_count, 
            encoding="utf-8"
        )
        
        # Create console handler
        console_handler = logging.StreamHandler(sys.stdout)

        # Set format
        fmt = fmt or "%(asctime)s %(levelname)s [%(name)s] %(message)s"
        formatter = logging.Formatter(fmt)
        file_handler.setFormatter(formatter)
        console_handler.setFormatter(formatter)

        # Add handlers to logger
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)

        cls._logger = logger
        cls._configured = True

    @classmethod
    def _get_logger(cls) -> logging.Logger:
        if not cls._configured or cls._logger is None:
            cls.configure()
        return cls._logger  # type: ignore[return-value]

    @classmethod
    def _safe_log(cls, log_func, msg: Any, *args, **kwargs) -> None:
        """Safe logging wrapper that handles various argument formats."""
        # Convert msg to string if it's not already
        if not isinstance(msg, str):
            msg = str(msg)
        
        # If we have args but msg doesn't contain format placeholders,
        # combine msg with args into a single message
        if args and not any(placeholder in msg for placeholder in ['%', '{', '}']):
            # Convert all args to strings and join them
            arg_strings = [str(arg) for arg in args]
            full_message = msg + ' ' + ' '.join(arg_strings)
            log_func(full_message, **kwargs)
        else:
            # Normal logging with format strings
            log_func(msg, *args, **kwargs)

    @classmethod
    def debug(cls, msg: Any, *args, **kwargs) -> None:
        cls._safe_log(cls._get_logger().debug, msg, *args, **kwargs)

    @classmethod
    def info(cls, msg: Any, *args, **kwargs) -> None:
        cls._safe_log(cls._get_logger().info, msg, *args, **kwargs)

    @classmethod
    def warning(cls, msg: Any, *args, **kwargs) -> None:
        cls._safe_log(cls._get_logger().warning, msg, *args, **kwargs)

    @classmethod
    def error(cls, msg: Any, *args, **kwargs) -> None:
        cls._safe_log(cls._get_logger().error, msg, *args, **kwargs)

    @classmethod
    def critical(cls, msg: Any, *args, **kwargs) -> None:
        cls._safe_log(cls._get_logger().critical, msg, *args, **kwargs)

    @classmethod
    def exception(cls, msg: Any, *args, exc_info: bool = True, **kwargs) -> None:
        cls._safe_log(cls._get_logger().exception, msg, *args, exc_info=exc_info, **kwargs)


# Module-level convenience wrappers for direct import
def debug(msg: Any, *args, **kwargs) -> None:
    LoggerHelper.debug(msg, *args, **kwargs)


def info(msg: Any, *args, **kwargs) -> None:
    LoggerHelper.info(msg, *args, **kwargs)


def warning(msg: Any, *args, **kwargs) -> None:
    LoggerHelper.warning(msg, *args, **kwargs)


def error(msg: Any, *args, **kwargs) -> None:
    LoggerHelper.error(msg, *args, **kwargs)


def critical(msg: Any, *args, **kwargs) -> None:
    LoggerHelper.critical(msg, *args, **kwargs)


def exception(msg: Any, *args, **kwargs) -> None:
    LoggerHelper.exception(msg, *args, **kwargs)


# Export LoggerHelper as Logger for backward compatibility
Logger = LoggerHelper

__all__ = ["Logger", "LoggerHelper", "debug", "info", "warning", "error", "critical", "exception"]