import asyncio
from contextlib import asynccontextmanager
from typing import Optional, AsyncGenerator

import aiomysql
from aiomysql import Pool, Connection
from app.configs.settings import settings
from app.utils.LoggerHelper import LoggerHelper

class MySQLPool:
    """异步 MySQL 数据库连接池管理（线程/协程安全单例）"""
    _pool: Optional[Pool] = None
    _lock = asyncio.Lock()
    
    @classmethod
    async def get_pool(cls) -> Pool:
        """获取数据库连接池（懒加载 + 并发安全）"""
        if cls._pool is None:
            async with cls._lock:
                if cls._pool is None:  # 双重检查锁
                    await cls._create_pool()
        return cls._pool

    @classmethod
    async def _create_pool(cls) -> None:
        """创建数据库连接池"""
        try:
            cls._pool = await aiomysql.create_pool(
                host=settings.DB_HOST,
                port=settings.DB_PORT,
                user=settings.DB_USER,
                password=settings.DB_PASSWORD,
                db=settings.DB_NAME,
                minsize=settings.DB_POOL_SIZE,
                maxsize=settings.DB_POOL_SIZE * 2,
                autocommit=False,  # 手动控制事务
                cursorclass=aiomysql.DictCursor,
                charset='utf8mb4',
                use_unicode=True,
                pool_recycle=settings.DB_POOL_RECYCLE,
                echo=False,
            )
            LoggerHelper.info("Database connection pool created successfully.")
        except Exception as e:
            LoggerHelper.error(f"Failed to create database connection pool: {e}")
            raise

    @classmethod
    @asynccontextmanager
    async def get_connection(cls) -> AsyncGenerator[Connection, None]:
        """从连接池中获取一个连接（自动释放）"""
        pool = await cls.get_pool()
        conn = await pool.acquire()
        try:
            yield conn
        finally:
            pool.release(conn)

    @classmethod
    async def close_pool(cls) -> None:
        """关闭连接池（应在应用关闭时调用）"""
        if cls._pool is not None:
            cls._pool.close()
            await cls._pool.wait_closed()
            cls._pool = None
            LoggerHelper.info("Database connection pool closed.")


# ========================
# FastAPI 依赖注入函数
# ========================

async def get_db_connection() -> AsyncGenerator[Connection, None]:
    """
    获取原始数据库连接（不自动提交事务）。
    适用于需要手动控制事务或只读操作的场景。
    """
    async with MySQLPool.get_connection() as conn:
        yield conn

async def get_db_transaction() -> AsyncGenerator[Connection, None]:
    """
    获取数据库连接并自动管理事务：
      - 正常结束 → 自动 commit
      - 发生异常 → 自动 rollback
    适用于写操作（如创建、更新、删除）。
    """
    async with MySQLPool.get_connection() as conn:
        try:
            yield conn
        except Exception:
            await conn.rollback()
            LoggerHelper.exception("Rollback due to exception in transaction.")
            raise
        else:
            await conn.commit()