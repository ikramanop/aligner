pub const INIT_QUERIES: &[&str] = &[
    "create table if not exists base_matrices (
        id int primary key auto_increment,
        dim smallint not null,
        matrix_json json not null unique
    )",
    "create table if not exists align_tasks (
        id int primary key auto_increment,
        hash text not null unique,
        query_sequence_id text not null,
        query_sequence text not null,
        target_sequence_id text not null,
        target_sequence text not null,
        kd_value double not null,
        r_squared_value double not null,
        del_value double not null,
        dim_value smallint not null,
        matrices_volume_value smallint not null,
        status text not null,
        p_value double
    )",
    "create table if not exists align_subtasks (
        id int primary key auto_increment,
        task_id int not null,
        f_value double not null,
        matrix_json json not null,
        result_query_sequence text,
        result_target_sequence text,
        foreign key (task_id) references align_tasks (id)
    )",
    "create table if not exists result_matrices (
        id int primary key auto_increment,
        task_id int not null,
        f_value double not null,
        matrix_json json not null,
        result_query_sequence text,
        result_target_sequence text,
        foreign key (task_id) references align_tasks (id)
    )",
];

pub const GET_BASE_MATRICES_WITH_LIMIT: &str = "
select matrix_json from base_matrices
where dim = ?
limit ?
";

pub const INSERT_BASE_MATRIX: &str = "
insert into base_matrices (dim, matrix_json) values (?, ?)
";

pub const INSERT_ALIGN_TASK: &str = "
insert into align_tasks (hash, query_sequence_id, query_sequence, target_sequence_id, target_sequence, kd_value, r_squared_value, del_value, dim_value, matrices_volume_value, status)
values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
";

pub const GET_ALIGN_TASK_ID_BY_HASH: &str = "
select id from align_tasks
where hash = ?
";

pub const INSERT_ALIGN_SUBTASKS: &str = "
insert into align_subtasks (task_id, f_value, matrix_json, result_query_sequence, result_target_sequence)
values (?, ?, ?, ?, ?)
";

pub const GET_PERCENTAGE_BY_HASH: &str = "
select count(1) / matrices_volume_value * 100 percentage from align_tasks at
    inner join align_subtasks a on at.id = a.task_id
group by hash
having hash = ?
";

pub const GET_SUBTASK_WITH_MAX_F_VALUE_BY_HASH: &str = "
select f_value, matrix_json, result_query_sequence, result_target_sequence from align_subtasks a
    inner join align_tasks at on a.task_id = at.id
    where hash = ?
order by f_value
limit 1
";

pub const INSERT_RESULT_MATRIX: &str = "
insert into result_matrices (task_id, f_value, matrix_json, result_query_sequence, result_target_sequence)
values (?, ?, ?, ?, ?)
";

pub const DELETE_SUBTASKS_BY_HASH: &str = "
delete a from align_subtasks a
    inner join align_tasks at on a.task_id = at.id
    where at.hash = ?
";

pub const GET_IDS_WITH_NULL_P_VALUE: &str = "
select id from align_tasks where p_value is null
";

pub const GET_RESULT_MATRIX_BY_TASK_ID: &str = "
select query_sequence, target_sequence, f_value, del_value, matrix_json from result_matrices rm
    inner join align_tasks at on rm.task_id = at.id
    where task_id = ?
";

pub const ADD_P_VALUE_BY_ID: &str = "
update align_tasks set p_value = ? where id = ?
";

pub const GET_ALL_HASHES: &str = "
select hash from align_tasks where p_value is null
";
