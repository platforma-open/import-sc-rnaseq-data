import type { InferOutputsType } from '@platforma-sdk/model';
import { BlockModel } from '@platforma-sdk/model';

export type BlockArgs = {
  name?: string;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({})

  .output('tengoMessage', (ctx) => ctx.outputs?.resolve('tengoMessage')?.getDataAsJson())

  .sections((_ctx) => [{ type: 'link', href: '/', label: 'Main' }])

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
