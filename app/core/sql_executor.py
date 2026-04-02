from typing import Optional, List, Any, Dict
from aiomysql import Connection
from app.utils.LoggerHelper import LoggerHelper

class SQLExecutor:
    """
    简单 SQL 执行器（仅支持无参数的原始 SQL 字符串）
    """

    def __init__(self, conn: Connection):
        self.conn = conn

    async def get_one(self, sql: str) -> Optional[Dict[str, Any]]:
        """执行查询，返回第一行结果（字典），无结果则返回 None"""
        LoggerHelper.debug(f"Executing get_one: {sql}")
        async with self.conn.cursor() as cursor:
            await cursor.execute(sql)
            row = await cursor.fetchone()
            return row  # aiomysql.DictCursor 返回 dict 或 None

    async def get_value(self, sql: str) -> Any:
        """执行查询，返回第一行第一列的值（如 COUNT(*)）"""
        LoggerHelper.debug(f"Executing get_value: {sql}")
        async with self.conn.cursor() as cursor:
            await cursor.execute(sql)
            row = await cursor.fetchone()
            if row is None:
                return None
            return next(iter(row.values()))  # 取第一个字段的值

    async def get_list(self, sql: str) -> List[Dict[str, Any]]:
        """执行查询，返回所有结果（列表 of 字典）"""
        LoggerHelper.debug(f"Executing get_list: {sql}")
        async with self.conn.cursor() as cursor:
            await cursor.execute(sql)
            rows = await cursor.fetchall()
            return rows or []

    async def execute_sql(self, sql: str) -> int:
        """执行非查询语句（INSERT/UPDATE/DELETE），返回受影响行数"""
        LoggerHelper.debug(f"Executing execute_sql: {sql}")
        async with self.conn.cursor() as cursor:
            await cursor.execute(sql)
            return cursor.rowcount