# Configure the AWS provider
provider "aws" {
  region = "eu-west-2"
  shared_credentials_files = ["$HOME/.aws/credentials"]
}

# Locals
# Defined in .env file in root directory. Must contain:
# - AWS_ACCOUNT_ID
# - AWS_REGION
# - ECR_REGISTRY
# - ECR_REPOSITORY
# - ECR_IMAGE_TAG
# - S3_BUCKET_NAME
# - S3_DATA_PREFIX
# - S3_EXECUTION_FOLDER
# - ACCESS_KEY_PARAMETER_NAME
# - SECRET_KEY_PARAMETER_NAME
locals {
  envs = { for tuple in regexall("(.*)=(.*)", file(".env")) : tuple[0] => sensitive(tuple[1]) }
}

## Network resources

# Create the VPC
resource "aws_vpc" "shrnaseq_vpc" {
  cidr_block = "10.0.0.0/16"

  tags = {
    Name = "shrnaseq VPC"
  }
}

resource "aws_security_group" "shrnaseq_security_group" {
  name        = "shrnaseq_security_group"
  description = "Security group that allows all outbound traffic"
  vpc_id      = aws_vpc.shrnaseq_vpc.id

  ingress {
   protocol         = "tcp"
   from_port        = 80
   to_port          = 80
   cidr_blocks      = ["0.0.0.0/0"]
   ipv6_cidr_blocks = ["::/0"]
  }
 
  ingress {
   protocol         = "tcp"
   from_port        = 443
   to_port          = 443
   cidr_blocks      = ["0.0.0.0/0"]
   ipv6_cidr_blocks = ["::/0"]
  }
 
  egress {
   protocol         = "-1"
   from_port        = 0
   to_port          = 0
   cidr_blocks      = ["0.0.0.0/0"]
   ipv6_cidr_blocks = ["::/0"]
  }
}

# Create the public subnet
resource "aws_subnet" "shrnaseq_public_subnet" {
  vpc_id            = aws_vpc.shrnaseq_vpc.id
  cidr_block        = "10.0.1.0/24"
  map_public_ip_on_launch = true

  tags = {
    Name = "shrnaseq Public Subnet"
  }
}

# Create the Internet gateway
resource "aws_internet_gateway" "shrnaseq_igw" {
  vpc_id = aws_vpc.shrnaseq_vpc.id

  tags = {
    Name = "shrnaseq Internet Gateway"
  }
}

# Create the route table for the public subnet
resource "aws_route_table" "shrnaseq_public_route_table" {
  vpc_id = aws_vpc.shrnaseq_vpc.id

  route {
    cidr_block = "0.0.0.0/0"
    gateway_id = aws_internet_gateway.shrnaseq_igw.id
  }

  tags = {
    Name = "shrnaseq Public Route Table"
  }
}

# Associate the public subnet with the public route table
resource "aws_route_table_association" "shrnaseq_public_subnet_association" {
  subnet_id      = aws_subnet.shrnaseq_public_subnet.id
  route_table_id = aws_route_table.shrnaseq_public_route_table.id
}

## IAM

# Create the policy which allows other actions for the EC2 instance
data "aws_iam_policy_document" "cloudwatch_domain_join_policy" {
  statement {
    actions = [
      "ecs:RunTask",
      "s3:ListAllMyBuckets",
      "s3:*"
    ]
    effect = "Allow"
    resources = ["*"]
  }
  statement {
    actions   = ["iam:PassRole"]
    effect    = "Allow"
    resources = ["*"]
    condition {
      test     = "StringLike"
      variable = "iam:PassedToService"
      values   = ["ecs-tasks.amazonaws.com"]
    }
  }
}

data "aws_iam_policy_document" "ecs_assume_role_policy" {
  statement {
    actions = [
      "sts:AssumeRole"
    ]
    effect = "Allow"
    principals {
      type        = "Service"
      identifiers = ["ecs-tasks.amazonaws.com"]
    }
  }
}

data "aws_iam_policy_document" "ecs_task_policy" {
  statement {
    actions = ["ssm:GetParameters", "secretsmanager:GetSecretValue"]
    effect = "Allow"
    resources = [
      "arn:aws:ssm:${local.envs["AWS_REGION"]}:${local.envs["AWS_ACCOUNT_ID"]}:parameter/${local.envs["ACCESS_KEY_PARAMETER_NAME"]}", 
      "arn:aws:ssm:${local.envs["AWS_REGION"]}:${local.envs["AWS_ACCOUNT_ID"]}:parameter/${local.envs["SECRET_KEY_PARAMETER_NAME"]}"
    ]
  }
  statement {
    actions = ["logs:GetLogEvents", "logs:PutLogEvents", "logs:CreateLogStream", "logs:DescribeLogStreams", "logs:PutRetentionPolicy", "logs:CreateLogGroup"]
    effect = "Allow"
    resources = ["arn:aws:logs:${local.envs["AWS_REGION"]}:${local.envs["AWS_ACCOUNT_ID"]}:log-group:awslogs-wordpress:log-stream:*"]
  }
}

# Create an IAM role for the ECS task execution
resource "aws_iam_role" "ecs_task_execution_role" {
  name = "ecs_task_execution_role"

  assume_role_policy = data.aws_iam_policy_document.ecs_assume_role_policy.json

  inline_policy {
    name = "ecs_task_role_policy"
    policy = data.aws_iam_policy_document.ecs_task_policy.json
  }
}

data "aws_iam_policy_document" "cloudwatch_assume_role_policy" {
  statement {
    actions = [
      "sts:AssumeRole"
    ]
    effect = "Allow"
    principals {
      type        = "Service"
      identifiers = ["events.amazonaws.com"]
    }

  }
}

# Create an IAM role for CloudWatch Events to assume
resource "aws_iam_role" "cloudwatch_role" {
  name = "cloudwatch_role"

  assume_role_policy = data.aws_iam_policy_document.cloudwatch_assume_role_policy.json

  inline_policy {
    name = "cloudwatch_role_policy"
    policy = data.aws_iam_policy_document.cloudwatch_domain_join_policy.json
  }
}

# Attach the CloudWatchEventsBuiltInTargetExecutionAccess to the cloudwatch IAM role
resource "aws_iam_policy_attachment" "cloudwatch_target_execution_role_policy" {
  name       = "CloudWatchEventsInvocationAccess"
  roles      = [aws_iam_role.cloudwatch_role.name]
  policy_arn = "arn:aws:iam::aws:policy/service-role/CloudWatchEventsBuiltInTargetExecutionAccess"
}

# Attach the CloudWatchEventsInvocationAccess to the cloudwatch IAM role
resource "aws_iam_policy_attachment" "cloudwatch_events_invocation_role_policy" {
  name       = "cloudwatch_events_invocation_role_policy"
  roles      = [aws_iam_role.cloudwatch_role.name]
  policy_arn = "arn:aws:iam::aws:policy/service-role/CloudWatchEventsInvocationAccess"
}

# Attach the AmazonECSTaskExecutionRolePolicy to the ecs IAM role
resource "aws_iam_policy_attachment" "ecs_task_execution_role_policy" {
  name       = "ecs_task_execution_role_policy"
  roles      = [aws_iam_role.ecs_task_execution_role.name]
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}

# Attach the AmazonECSTaskExecutionRolePolicy to the cloudwatch IAM role
resource "aws_iam_policy_attachment" "cloudwatch_task_execution_role_policy" {
  name       = "cloudwatch_task_execution_role_policy"
  roles      = [aws_iam_role.cloudwatch_role.name]
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonECSTaskExecutionRolePolicy"
}

## Task definition

# Create the task definition
resource "aws_ecs_task_definition" "shrnaseq_task_definition" {
  family                   = "shrnaseq-task-definition"
  network_mode             = "awsvpc"
  requires_compatibilities = ["FARGATE"]
  cpu                      = var.ecs_task_cpu
  memory                  = var.ecs_task_memory
  execution_role_arn       = aws_iam_role.ecs_task_execution_role.arn
  task_role_arn            = aws_iam_role.ecs_task_execution_role.arn

  container_definitions = <<EOF
  [
    {
      "name": "snakemake-container",
      "image": "${local.envs["ECR_REGISTRY"]}/${local.envs["ECR_REPOSITORY"]}:${local.envs["ECR_IMAGE_TAG"]}",
      "essential": true,
      "memory": ${var.ecs_task_memory},
      "secrets": [
        {
          "name": "AWS_ACCESS_KEY_ID",
          "valueFrom": "arn:aws:ssm:${local.envs["AWS_REGION"]}:${local.envs["AWS_ACCOUNT_ID"]}:parameter/${local.envs["ACCESS_KEY_PARAMETER_NAME"]}"
        },
        {
          "name": "AWS_SECRET_ACCESS_KEY",
          "valueFrom": "arn:aws:ssm:${local.envs["AWS_REGION"]}:${local.envs["AWS_ACCOUNT_ID"]}:parameter/${local.envs["SECRET_KEY_PARAMETER_NAME"]}"
        }
      ],
      "command": [
        "snakemake", 
        "--use-conda", 
        "--cores", 
        "all", 
        "--default-remote-provider", 
        "S3", 
        "--default-remote-prefix", 
        "${local.envs["S3_BUCKET_NAME"]}/${local.envs["S3_DATA_PREFIX"]}", 
        "--envvars", 
        "AWS_ACCESS_KEY_ID", 
        "AWS_SECRET_ACCESS_KEY"
      ],
      "logConfiguration": {
        "logDriver": "awslogs",
        "options": {
          "awslogs-create-group": "true",
          "awslogs-group": "awslogs-wordpress",
          "awslogs-region": "${local.envs["AWS_REGION"]}",
          "awslogs-stream-prefix": "awslogs-shrnaseq"
        }
      }
    }
  ]
  EOF
}

## Event rule

# Create the event rule
resource "aws_cloudwatch_event_rule" "shrnaseq_event_rule" {
  name        = "shrnaseq-event-rule"
  description = "Event rule that runs a Fargate task when a file is uploaded to an S3 bucket"
  event_pattern = <<EOF
  { "source": [ "aws.s3" ], 
    "detail-type": [ "Object Created" ], 
    "detail": { 
      "bucket": { 
        "name": [ 
          "${local.envs["S3_BUCKET_NAME"]}"
        ] 
      },
      "object": { 
        "key": [{"prefix" : "${local.envs["S3_DATA_PREFIX"]}/${local.envs["S3_EXECUTION_FOLDER"]}"}] 
      }
    } 
  }
  EOF
}

# Create the ECS target for the event rule
resource "aws_cloudwatch_event_target" "shrnaseq_event_target" {
  rule      = aws_cloudwatch_event_rule.shrnaseq_event_rule.name
  arn       = aws_ecs_cluster.shrnaseq_cluster.arn
  role_arn  = aws_iam_role.cloudwatch_role.arn
  ecs_target {
    task_definition_arn = aws_ecs_task_definition.shrnaseq_task_definition.arn
    launch_type = "FARGATE"
    network_configuration {
      subnets         = [aws_subnet.shrnaseq_public_subnet.id]
      security_groups = [aws_security_group.shrnaseq_security_group.id]
      assign_public_ip = true
    }
  }
}

## ECS cluster
resource "aws_ecs_cluster" "shrnaseq_cluster" {
  name = "shrnaseq-cluster"

  tags = {
    Name = "shrnaseq ECS Cluster"
  }
}