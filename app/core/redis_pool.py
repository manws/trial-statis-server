from typing import Optional
import asyncio
import logging

import redis.asyncio as redis
from app.configs.settings import settings
from app.utils.LoggerHelper import LoggerHelper

class RedisPool:
    """Enhanced async Redis client holder with password support and health checks.
    
    Features:
    - Password authentication support
    - Connection timeout configuration
    - Health check mechanism
    - Auto-reconnect capability
    - Thread-safe singleton pattern

    Example:
        # Basic usage
        client = await RedisPool.get_client()
        await client.set('key', 'value')
        
        # Health check
        is_healthy = await RedisPool.ping()
        
        # Reconnect if needed
        await RedisPool.reconnect()
    """
    _client: Optional[redis.Redis] = None
    _lock = asyncio.Lock()

    @classmethod
    async def init_pool(cls) -> None:
        """初始化Redis连接池"""
        if cls._client is not None:
            return
            
        async with cls._lock:
            if cls._client is not None:
                return
                
            try:
                redis_url = f"redis://:{settings.REDIS_PASSWORD}@{settings.REDIS_HOST}:{settings.REDIS_PORT}/{settings.REDIS_DB}"
                cls._client = redis.from_url(
                    redis_url,
                    decode_responses=True,
                    socket_connect_timeout=settings.REDIS_SOCKET_CONNECT_TIMEOUT,
                    socket_timeout=settings.REDIS_SOCKET_TIMEOUT,
                    retry_on_timeout=settings.REDIS_RETRY_ON_TIMEOUT,
                    health_check_interval=settings.REDIS_HEALTH_CHECK_INTERVAL
                )
                
                # 验证连接
                await cls._client.ping()
            except Exception as e:
                LoggerHelper.error(f"初始化Redis连接失败: {e}")
                cls._client = None
                raise

    @classmethod
    async def get_client(cls) -> redis.Redis:
        """获取Redis客户端实例"""
        if cls._client is None:
            await cls.init_pool()
        return cls._client  # type: ignore[return-value]

    @classmethod
    async def ping(cls) -> bool:
        """检查Redis连接健康状态"""
        try:
            client = await cls.get_client()
            result = await client.ping()
            return result
        except Exception as e:
            LoggerHelper.warning(f"Redis health check failed: {e}")
            return False

    @classmethod
    async def reconnect(cls) -> None:
        """重新连接Redis"""
        LoggerHelper.info("Attempting to reconnect Redis...")
        await cls.close_pool()
        await cls.init_pool()
        LoggerHelper.info("Redis reconnection completed")

    @classmethod
    async def close_pool(cls) -> None:
        """关闭Redis连接池"""
        if cls._client is None:
            return
            
        try:
            LoggerHelper.info("Closing Redis connection pool...")
            # 关闭客户端连接
            await cls._client.close()
        except Exception as e:
            LoggerHelper.error(f"Error closing Redis client: {e}")
            # 最佳努力断开连接
            try:
                await cls._client.connection_pool.disconnect()
            except Exception as disconnect_error:
                LoggerHelper.error(f"Error disconnecting Redis pool: {disconnect_error}")
        finally:
            cls._client = None
            LoggerHelper.info("Redis connection pool closed")


__all__ = ["RedisPool"]